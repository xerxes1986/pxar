#include "datapipe.h"
#include "helper.h"
#include "log.h"
#include "constants.h"
#include "exceptions.h"

namespace pxar {

  rawEvent* dtbEventSplitter::SplitDeser400() {
    record.Clear();

    // If last one had Event end marker, get a new sample:
    if (!nextStartDetected) { Get(); }

    // If new sample does not have start marker keep on reading until we find it:
    if ((GetLast() & 0xe000) != 0xa000) {
      record.SetStartError();
      Get();
    }
    record.Add(GetLast());

    // Else keep reading and adding samples until we find any marker.
    while ((Get() & 0xe000) != 0xc000) {
      // Check if the last read sample has Event end marker:
      if ((GetLast() & 0xe000) == 0xa000) {
	record.SetEndError();
	nextStartDetected = true;
	return &record;
      }
      // If total Event size is too big, break:
      if (record.GetSize() < 40000) record.Add(GetLast());
      else record.SetOverflow();
    }
    record.Add(GetLast());
    nextStartDetected = false;

    LOG(logDEBUGPIPES) << "SINGLE SPLIT EVENT:";
    LOG(logDEBUGPIPES) << listVector(record.data,true);
    LOG(logDEBUGPIPES) << "-------------------------";

    return &record;
  }

  rawEvent* dtbEventSplitter::SplitDeser160() {
    record.Clear();

    // If last one had Event end marker, get a new sample:
    if (GetLast() & 0x4000) { Get(); }

    // If new sample does not have start marker keep on reading until we find it:
    if (!(GetLast() & 0x8000)) {
      record.SetStartError();
      while (!(GetLast() & 0x8000)) Get();
    }

    // Else keep reading and adding samples until we find any marker.
    do {
      // If total Event size is too big, break:
      if (record.GetSize() >= 40000) {
	record.SetOverflow();
	break;
      }

      // FIXME Very first Event starts with 0xC - which srews up empty Event detection here!
      // If the Event start sample is also Event end sample, write and quit:
      if((GetLast() & 0xc000) == 0xc000) { break; }

      record.Add(GetLast());
    } while ((Get() & 0xc000) == 0);

    // Check if the last read sample has Event end marker:
    if (GetLast() & 0x4000) record.Add(GetLast());
    // Else set Event end error:
    else record.SetEndError();

    LOG(logDEBUGPIPES) << "SINGLE SPLIT EVENT:";
    if(GetDeviceType() < ROC_PSI46DIG) { LOG(logDEBUGPIPES) << listVector(record.data,false,true); }
    else { LOG(logDEBUGPIPES) << listVector(record.data,true); }
    LOG(logDEBUGPIPES) << "-------------------------";

    return &record;
  }

  rawEvent* dtbEventSplitter::SplitSoftTBM() {
   record.Clear();

   // If last one had Event end marker, get a new sample:
   if (!nextStartDetected) { Get(); }

   // If new sample does not have start marker keep on reading until we find it:
   while ((GetLast() & 0xe000) != 0xa000) {
     record.SetStartError();
     Get();
   }
   Get();
   //record.Add(GetLast());

   // Else keep reading and adding samples until we find any trailer marker.
   while ((Get() & 0xe000) != 0xe000) {
     // Check if the last read sample has Event end marker:
     if ((GetLast() & 0xe000) == 0xa000) {
       record.SetEndError();
       nextStartDetected = true;
       return &record;
     }

     // If total Event size is too big, break:
     if (record.GetSize() < 40000) record.Add(GetLast());
     else record.SetOverflow();
   }

   LOG(logDEBUGPIPES) << "-------------------------";
   LOG(logDEBUGPIPES) << listVector(record.data,true);

   return &record;
  }

  rawEvent* passthroughSplitter::Read() {
    record.Clear();
    try {
      do {
	record.Add(Get());
      } while(1);
    }
    catch(dsBufferEmpty) {}
    return &record;
  }

  void dtbEventDecoder::CheckInvalidWord(uint16_t v) {
    // Check last bit of identifier nibble to be zero:
    if((v & 0x1000) == 0x0000) { return; }
    decodingStats.m_errors_event_invalid_words++;
  }

  void dtbEventDecoder::CheckEventID(uint16_t v) {
    // After startup, register the first event ID:
    if(eventID == -1) { eventID = (v&0x00ff); }

    // Check if TBM event ID matches with expectation:
    if((v&0x00ff) != (eventID%256)) {
      LOG(logERROR) << "   Event ID mismatch:  local ID (" << static_cast<int>(eventID) 
		    << ") !=  TBM ID (" << static_cast<int>(v&0x00ff) << ")";
      decodingStats.m_errors_tbm_eventid_mismatch++;
      // To continue readout, set event ID to the currently decoded one:
      eventID = (v&0x00ff);
    }

    // Increment event counter:
    eventID = (eventID%256) + 1;
  }

  Event* dtbEventDecoder::DecodeDeser400() {

    roc_Event.Clear();
    rawEvent *sample = Get();

    // Count possibe error states:
    if(sample->IsStartError()) { decodingStats.m_errors_event_start++; }
    if(sample->IsEndError()) { decodingStats.m_errors_event_stop++; }
    if(sample->IsOverflow()) { decodingStats.m_errors_event_overflow++; }

    unsigned int raw = 0;
    unsigned int pos = 0;
    unsigned int size = sample->GetSize();
    decodingStats.m_info_words_read += size;

    uint16_t v;
    bool tmpError = false;

    // Count the ROC headers:
    int16_t roc_n = -1;

    // Check if ROC has inverted pixel address (ROC_PSI46DIG):
    bool invertedAddress = ( GetDeviceType() == ROC_PSI46DIG ? true : false );
    

    // --- decode TBM header ---------------------------------

    // TBM Header 1:
    v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
    CheckInvalidWord(v);
    if ((v & 0xe000) != 0xa000) tmpError = true;
    raw = (v & 0x00ff) << 8;
    LOG(logDEBUGPIPES) << "TBM " << static_cast<int>(GetChannel()) << " Evt ID " << static_cast<int>(v&0x00ff);

    // Check for correct TBM event ID:
    CheckEventID(v);

    // TBM Header 2:
    v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
    CheckInvalidWord(v);
    if ((v & 0xe000) != 0x8000) tmpError = true;
    raw += v & 0x00ff;
    LOG(logDEBUGPIPES) << "\t Data ID " << static_cast<int>(((v & 0x00c0) >> 6)) 
		       << " Value " << static_cast<int>((v & 0x003f));

    if(tmpError) { decodingStats.m_errors_tbm_header++; }
    tmpError = false;

    roc_Event.header = raw;

    // --- decode ROC data -----------------------------------

    // while ROC header
    v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
    CheckInvalidWord(v);

    while ((v & 0xe000) == 0x4000) { // ROC Header

      // Count ROC Headers up:
      roc_n++;

      // Check for DESER400 failure:
      if((v&0x0ff0) == 0x0ff0) {
	LOG(logCRITICAL) << "TBM " << static_cast<int>(GetChannel())
			 << " ROC " << static_cast<int>(roc_n)
			 << " header reports DESER400 failure!";
	decodingStats.m_errors_event_invalid_xor++;
	throw DataDecodingError("Invalid XOR eye diagram encountered.");
      }

      // Decode the readback bits in the ROC header:
      if(GetDeviceType() >= ROC_PSI46DIGV2) { evalReadback(static_cast<uint8_t>(roc_n),v); }

      v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
      CheckInvalidWord(v);
      while ((v & 0xe000) <= 0x2000) { // R0 ... R1

	for (int i = 0; i <= 1; i++) {
	  if ((v >> 13) != i) { // R<i>
	    if (v & 0x8000) { // TBM header/trailer
	      // Unexpected arrival of TBM marker - pixel data is incomplete:
	      decodingStats.m_errors_pixel_incomplete++;
	      v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
	      CheckInvalidWord(v);
	      goto trailer;
	    }
	  }
	  raw = (raw << 12) + (v & 0x0fff);
	  v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
	  CheckInvalidWord(v);
	}

	try {
	  // Check if this is just fill bits of the TBM09 data stream 
	  // accounting for the other channel:
	  if(GetTokenChainLength() == 4 && (raw&0xffffff) == 0xffffff) {
	    LOG(logDEBUGPIPES) << "Empty hit detected (TBM09 data streams). Skipping.";
	    continue;
	  }

	  // Get the correct ROC id: Channel number x ROC offset (= token chain length)
	  // TBM08x: channel 0: 0-7, channel 1: 8-15
	  // TBM09x: channel 0: 0-3, channel 1: 4-7, channel 2: 8-11, channel 3: 12-15
	  pixel pix(raw,static_cast<uint8_t>(roc_n + GetChannel()*GetTokenChainLength()),invertedAddress);
	  roc_Event.pixels.push_back(pix);
	  decodingStats.m_info_pixels_valid++;
	}
	catch(DataInvalidAddressError /*&e*/){
	  // decoding of raw address lead to invalid address
	  decodingStats.m_errors_pixel_address++;
	}
	catch(DataInvalidPulseheightError /*&e*/){
	  // decoding of pulse height featured non-zero fill bit
	  decodingStats.m_errors_pixel_pulseheight++;
	}
	catch(DataCorruptBufferError /*&e*/){
	  // decoding returned row 80 - corrupt data buffer
	  decodingStats.m_errors_pixel_buffer_corrupt++;
	}
      }
      //if (roc.error) x.error |= 0x0001;
      //x.roc.push_back(roc);
    }

    // --- decode TBM trailer --------------------------------
  trailer:
    raw = 0;

    // T1
    if ((v & 0xe000) != 0xe000) tmpError = true;
    raw = (v & 0x00ff) << 8;
    LOG(logDEBUGPIPES) << "TBM " << static_cast<int>(GetChannel()) << " trailer content: " << std::hex << (v&0x00ff) << std::dec;
    LOG(logDEBUGPIPES) << "\t Token Pass " << ((v&0x0080) == 0x0080 ? "FALSE" : "TRUE");
    LOG(logDEBUGPIPES) << "\t REST " << static_cast<int>(v&0x0040);
    LOG(logDEBUGPIPES) << "\t RESR " << static_cast<int>(v&0x0020);
    LOG(logDEBUGPIPES) << "\t Sync Err " << static_cast<int>(v&00010);
    LOG(logDEBUGPIPES) << "\t Sync Trigger " << static_cast<int>(v&0x0008);
    LOG(logDEBUGPIPES) << "\t Clear Trig Count " << static_cast<int>(v&0x0004);
    LOG(logDEBUGPIPES) << "\t Cal Trigger " << static_cast<int>(v&0x0002);
    LOG(logDEBUGPIPES) << "\t Stack Full " << static_cast<int>(v&0x0001);

    // T2
    v = (pos < size) ? (*sample)[pos++] : 0x6000; //MDD_ERROR_MARKER;
    CheckInvalidWord(v);
    if ((v & 0xe000) != 0xc000) tmpError = true;
    raw += v & 0x00ff;
    LOG(logDEBUGPIPES) << "\t Stack Full Now " << static_cast<int>(v&0x0080);
    LOG(logDEBUGPIPES) << "\t PKAM Reset " << static_cast<int>(v&0x0040);
    LOG(logDEBUGPIPES) << "\t Stack Count " << static_cast<int>(v&0x003f);

    if(tmpError) { decodingStats.m_errors_tbm_trailer++; }

    roc_Event.trailer = raw;

    // Check event validity (empty, missing ROCs...):
    CheckEventValidity(roc_n);

    LOG(logDEBUGPIPES) << roc_Event;
    return &roc_Event;
  }

  void dtbEventDecoder::AverageAnalogLevel(int32_t &variable, int16_t dataword) {
    // Check if this is the initial measurement:
    if(variable > 0xff) { variable = expandSign(dataword & 0x0fff); }
    // Average the variable:
    else { variable = (variable + expandSign(dataword & 0x0fff))/2; }
  }

  Event* dtbEventDecoder::DecodeAnalog() {

    roc_Event.Clear();
    rawEvent *sample = Get();
    int16_t roc_n = -1;

    // Count possibe error states:
    if(sample->IsStartError()) { decodingStats.m_errors_event_start++; }
    if(sample->IsEndError()) { decodingStats.m_errors_event_stop++; }
    if(sample->IsOverflow()) { decodingStats.m_errors_event_overflow++; }

    unsigned int n = sample->GetSize();
    decodingStats.m_info_words_read += n;

    if (n >= 3) {
      // Reserve expected number of pixels from data length (subtract ROC headers):
      if (n - 3*GetTokenChainLength() > 0) roc_Event.pixels.reserve((n - 3*GetTokenChainLength())/6);

      unsigned int pos = 0;
      while (pos+3 <= n) {
	// Check if we have another ROC header (UB and B levels):
	// Here we have to assume the first two words are a ROC header because we rely on
	// its Ultrablack and Black level as initial values for auto-calibration:
	int16_t levelS = (black - ultrablack)/8;

	if(roc_n < 0 || 
	   // Ultrablack level:
	   ((ultrablack-levelS < expandSign((*sample)[pos] & 0x0fff) && ultrablack+levelS > expandSign((*sample)[pos] & 0x0fff))
	   // Black level:
	    && (black-levelS < expandSign((*sample)[pos+1] & 0x0fff) && black+levelS > expandSign((*sample)[pos+1] & 0x0fff)))) {

	  roc_n++;
	  // Save the lastDAC value:
	  evalLastDAC(roc_n, (*sample)[pos+2] & 0x0fff);
	  roc_Event.header = (*sample)[pos+2] & 0x0fff;

	  // Iterate to improve ultrablack and black measurement:
	  AverageAnalogLevel(ultrablack, (*sample)[pos] & 0x0fff);
	  AverageAnalogLevel(black, (*sample)[pos+1] & 0x0fff);

	  LOG(logDEBUGPIPES) << "ROC Header: "
			     << expandSign((*sample)[pos] & 0x0fff) << " (avg. " << ultrablack << ") (UB) "
			     << expandSign((*sample)[pos+1] & 0x0fff) << " (avg. " << black << ") (B) "
			     << expandSign((*sample)[pos+2] & 0x0fff) << " (lastDAC) ";
	  pos += 3;
	}
	// We have a pixel hit:
	else {
	  // Not enough data for a new pixel hit (six words):
	  if(pos > n-6) break;

	  std::vector<uint16_t> data;
	  data.push_back((*sample)[pos] & 0x0fff);
	  data.push_back((*sample)[pos+1] & 0x0fff);
	  data.push_back((*sample)[pos+2] & 0x0fff);
	  data.push_back((*sample)[pos+3] & 0x0fff);
	  data.push_back((*sample)[pos+4] & 0x0fff);
	  data.push_back((*sample)[pos+5] & 0x0fff);
 
	  try{
	    LOG(logDEBUGPIPES) << "Trying to decode pixel: " << listVector(data,false,true);
	    pixel pix(data,roc_n,ultrablack,black);
	    roc_Event.pixels.push_back(pix);
	    decodingStats.m_info_pixels_valid++;
	  }
	  catch(DataDecodingError /*&e*/){
	    // decoding of raw address lead to invalid address
	    decodingStats.m_errors_pixel_address++;
	  }
	  // Advance read pointer by one pixel:
	  pos += 6;
	}
      }
    }

    // Check if empty or missing ROCs:
    CheckEventValidity(roc_n);

    LOG(logDEBUGPIPES) << roc_Event;
    return &roc_Event;
  }
  void dtbEventDecoder::CheckEventValidity(int16_t roc_n) {

    // Check that we found all expected ROC headers:
    // If the number of ROCs does not correspond to what we expect
    // clear the event and return:
    if(roc_n+1 != GetTokenChainLength()) {
      LOG(logERROR) << "Number of ROCs (" << static_cast<int>(roc_n+1)
		    << ") != Token Chain Length (" << static_cast<int>(GetTokenChainLength()) << ")";
      decodingStats.m_errors_roc_missing++;
      // Clearing event content:
      roc_Event.Clear();
    }
    // Count empty events
    else if(roc_Event.pixels.empty()) { decodingStats.m_info_events_empty++; }
    // Count valid events
    else { decodingStats.m_info_events_valid++; }
  }

  Event* dtbEventDecoder::DecodeDeser160() {

    roc_Event.Clear();

    // Count the ROC headers:
    int16_t roc_n = -1;

    // Check if ROC has inverted pixel address (ROC_PSI46DIG):
    bool invertedAddress = ( GetDeviceType() == ROC_PSI46DIG ? true : false );

    rawEvent *sample = Get();

    // Count possibe error states:
    if(sample->IsStartError()) { decodingStats.m_errors_event_start++; }
    if(sample->IsEndError()) { decodingStats.m_errors_event_stop++; }
    if(sample->IsOverflow()) { decodingStats.m_errors_event_overflow++; }

    unsigned int n = sample->GetSize();
    decodingStats.m_info_words_read += n;

    if (n > 0) {
      // Reserve expected number of pixels from data length (subtract ROC headers):
      if (n-GetTokenChainLength() > 0) roc_Event.pixels.reserve((n-GetTokenChainLength())/2);

      unsigned int pos = 0;

      // Loop over the full data:
      while (pos < n) {
	// Check if we have a ROC header:
	if(((*sample)[pos] & 0x0ffc) == 0x07f8) {
	  roc_Event.header = (*sample)[pos] & 0x0fff;
	  roc_n++;
	  pos++;

	  // Decode the readback bits in the ROC header:
	  if(GetDeviceType() >= ROC_PSI46DIGV2) { evalReadback(roc_n,roc_Event.header); }
	}
	// We have a pixel:
	// Require that we found at least one ROC header:
	else if(roc_n >= 0) {
	  // Not enough data for a new pixel hit (two words):
	  if(pos >= n-1) break;

	  uint32_t raw = ((*sample)[pos++] & 0x0fff) << 12;
	  raw += (*sample)[pos++] & 0x0fff;
	  try {
	    pixel pix(raw,roc_n,invertedAddress);
	    roc_Event.pixels.push_back(pix);
	    decodingStats.m_info_pixels_valid++;
	  }
	  catch(DataInvalidAddressError /*&e*/) {
	    // decoding of raw address lead to invalid address
	    decodingStats.m_errors_pixel_address++;
	  }
	  catch(DataInvalidPulseheightError /*&e*/) {
	    // decoding of pulse height featured non-zero fill bit
	    decodingStats.m_errors_pixel_pulseheight++;
	  }
	  catch(DataCorruptBufferError /*&e*/) {
	    // decoding returned row 80 - corrupt data buffer
	    decodingStats.m_errors_pixel_buffer_corrupt++;
	  }
	}
	// No ROC header present, try with next word:
	else {
	  pos++;
	}
      }
    }

    // Check if empty or missing ROCs:
    CheckEventValidity(roc_n);

    LOG(logDEBUGPIPES) << roc_Event;
    return &roc_Event;
  }

  void dtbEventDecoder::evalLastDAC(uint8_t roc, uint16_t val) {
    // Check if we have seen this ROC already:
    if(readback.size() <= roc) readback.resize(roc+1);
    readback.at(roc).push_back(val);

    LOG(logDEBUGPIPES) << "Readback ROC "
		       << static_cast<int>(roc+GetChannel()*GetTokenChainLength())
		       << " " << static_cast<int>(expandSign(val & 0x0fff));
  }

  void dtbEventDecoder::evalReadback(uint8_t roc, uint16_t val) {
    // Check if we have seen this ROC already:
    if(shiftReg.size() <= roc) shiftReg.resize(roc+1,0);
    shiftReg.at(roc) <<= 1;
    if(val&1) shiftReg.at(roc)++;

    // Count this bit:
    if(count.size() <= roc) count.resize(roc+1,0);
    count.at(roc)++;

    if(val&2) { // start marker detected
      if (count.at(roc) == 16) {
	// Write out the collected data:
	if(readback.size() <= roc) readback.resize(roc+1);
	readback.at(roc).push_back(shiftReg.at(roc));

	LOG(logDEBUGPIPES) << "Readback ROC "
			   << static_cast<int>(roc+GetChannel()*GetTokenChainLength())
			   << " " << ((readback.at(roc).back()>>8)&0x00ff)
			   << " (0x" << std::hex << ((readback.at(roc).back()>>8)&0x00ff)
			   << std::dec << "): " << (readback.at(roc).back()&0xff)
			   << " (0x" << std::hex << (readback.at(roc).back()&0xff)
			   << std::dec << ")";
      }
      else {
	// If this is the first readback cycle of the ROC, ignore the mismatch:
	if(readback.size() <= roc || readback.at(roc).empty()) {
	  LOG(logDEBUGAPI) << "ROC " << static_cast<int>(roc)
			   << ": first readback marker after "
			   << count.at(roc) << " readouts. Ignoring error condition.";
	}
	else {
	  LOG(logWARNING) << "ROC " << static_cast<int>(roc)
			  << ": Readback start marker after "
			  << count.at(roc) << " readouts!";
	  decodingStats.m_errors_roc_readback++;
	}
      }
      // Reset the counter for this ROC:
      count.at(roc) = 0;
    }
  }

  statistics dtbEventDecoder::getStatistics() { 
    // Automatically clear the statistics after it was read out:
    statistics tmp = decodingStats;
    decodingStats.clear();
    return tmp;
  }

  std::vector<std::vector<uint16_t> > dtbEventDecoder::getReadback() {
    // Automatically clear the readback vector after it was read out:
    std::vector<std::vector<uint16_t> > tmp = readback;
    readback.clear();
    return tmp;
  }
}
