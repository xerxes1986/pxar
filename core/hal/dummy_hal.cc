#include "hal.h"
#include "log.h"
#include "timer.h"
#include "helper.h"
#include "constants.h"
#include <fstream>
#include <stdlib.h>

using namespace pxar;

pxar::pixel getNoiseHit(uint8_t rocid, size_t i, size_t j) {

  // Generate a slightly random pulse height between 80 and 100:
  uint16_t pulseheight = rand()%20 + 80;

  // We can't pulse the same pixel twice in one trigger:
  size_t col = rand()%52;
  while(col == i) col = rand()%52;
  size_t row = rand()%80;
  while(row == j) row = rand()%80;

  pixel px = pixel(rocid,col,row,pulseheight);
  LOG(logDEBUGPIPES) << "Adding noise hit: " << px;
  return px;
}

bool isInTornadoRegion(size_t dac1min, size_t dac1max, size_t dac1, size_t dac2min, size_t dac2max, size_t dac2) {

  size_t epsilon = 5;
  double tornadowidth = 40;
  double ymax = dac2max / (1 + exp(0.05*((dac1max-dac1min+1)/2 - tornadowidth - dac1))) + dac2min;
  double ymin = dac2max / (1 + exp(0.05*((dac1max-dac1min+1)/2 + tornadowidth - dac1))) + dac2min;
  if(dac2 < ymax && dac2 > ymin) {
    size_t dymax = ymax - dac2;
    size_t dymin = dac2 - ymin;
    if(dymax < epsilon) return (rand()%(epsilon-dymax) == 0);
    else if(dymin < epsilon) return (rand()%(epsilon-dymin) == 0);
    else return true;
  }
  else return false;
}

void fillEvent(pxar::Event * evt, uint8_t rocid, size_t col, size_t row, uint32_t flags) {

  // Generate a slightly random pulse height between 90 and 100:
  uint16_t pulseheight = rand() % 2 + 90;

  // Introduce some address encoding issues:
  if((flags&FLAG_CHECK_ORDER) != 0 && col == 0 && row == 1) { evt->pixels.push_back(pixel(rocid,col,row+1,pulseheight));} // PX 0,1 answers as PX 0,2
  else if((flags&FLAG_CHECK_ORDER) != 0 && col == 0 && row == 2) { } // PX 0,2 is dead
  else if((flags&FLAG_CHECK_ORDER) != 0 && col == 0 && row == 6) { evt->pixels.push_back(pixel(rocid,col,row+1,pulseheight));} // PX 0,6 answers as PX 0,7
  else { evt->pixels.push_back(pixel(rocid,col,row,pulseheight)); }

  // If the full chip is unmasked, add some noise hits:
  if((flags&FLAG_FORCE_UNMASKED) != 0 && (rand()%4) == 0) { evt->pixels.push_back(getNoiseHit(rocid,col,row)); }

}

hal::hal(std::string /*name*/) :
  _initialized(false),
  _compatible(false),
  m_tbmtype(TBM_NONE),
  deser160phase(4),
  _currentTrgSrc(TRG_SEL_PG_DIR)
{
  // Print the useful SW/FW versioning info:
  PrintInfo();

  // Initialize rand():
  srand (time(NULL));

  // Check if all RPC calls are matched:
  if(CheckCompatibility()) {
    // Set compatibility flag
    _compatible = true;
  }
}

hal::~hal() {
}

bool hal::status() {
  return _initialized;
}

uint32_t hal::GetHashForString(const char * s)
{
  // Using some primes
  uint32_t h = 31;
  while (*s) { h = (h * 54059) ^ (s[0] * 76963); s++; }
  return h%86969;
}

uint32_t hal::GetHashForStringVector(const std::vector<std::string> & v)
{
  uint32_t ret = 0;
  for (size_t i=0; i < v.size(); i++) {ret += ((i+1)*(GetHashForString(v[i].c_str())));}
  return ret;
}

void hal::initTestboard(std::map<uint8_t,uint8_t> /*sig_delays*/, std::vector<std::pair<uint16_t,uint8_t> > /*pg_setup*/, uint16_t /*delaysum*/, double /*va*/, double /*vd*/, double /*ia*/, double /*id*/) {
  
  // We are ready for operations now, mark the HAL as initialized:
  _initialized = true;
}

void hal::SetupPatternGenerator(std::vector<std::pair<uint16_t,uint8_t> > /*pg_setup*/, uint16_t /*delaysum*/) {
}

void hal::setTestboardDelays(std::map<unsigned char, unsigned char, std::less<unsigned char>, std::allocator<std::pair<unsigned char const, unsigned char> > >) {
}

void hal::setTestboardPower(double, double, double, double) {
}

bool hal::flashTestboard(std::ifstream& /*flashFile*/) {
  LOG(logCRITICAL) << "ERROR UPGRADE: Could not upgrade this DTB version! It's a dummy after all!";
  return false;
}

void hal::initTBMCore(uint8_t type, std::map< uint8_t,uint8_t > /*regVector*/) {
  m_tbmtype = type;
}

void hal::setTBMType(uint8_t type) {
  LOG(logDEBUGHAL) << "Updating TBM type to 0x" << std::hex << type << std::dec << ".";
  m_tbmtype = type;
}

void hal::initROC(uint8_t roci2c, uint8_t type, std::map< uint8_t,uint8_t > dacVector) {
  // Set the pixel address inverted flag if we have the PSI46digV1 chip
  if(type == ROC_PSI46DIG || type == ROC_PSI46DIG_TRIG) {
    LOG(logDEBUGHAL) << "Pixel address is inverted in this ROC type.";
  }
  // FIXME
  m_roctype = type;

  // Programm all DAC registers according to the configuration data:
  LOG(logDEBUGHAL) << "Setting DAC vector for ROC@I2C " << static_cast<int>(roci2c) << ".";
  rocSetDACs(roci2c,dacVector);

}

void hal::PrintInfo() {
  LOG(logINFO) << "DTB startup information" << std::endl 
	       << "--- DTB info------------------------------------------" << std::endl
	       << " DUMMY DUMMY DUMMY " << std::endl
	       << "------------------------------------------------------";
}

bool hal::CheckCompatibility(){
  // We are though all checks, testboard is successfully connected:
  return true;
}

bool hal::FindDTB(std::string &/*usbId*/) {
  return true;
}

double hal::getTBia() {
  // Return the VA analog current in A:
  return (5.23);
}

double hal::getTBva(){
  // Return the VA analog voltage in V:
  return (5.23);
}

double hal::getTBid() {
  // Return the VD digital current in A:
  return (5.23);
}

double hal::getTBvd() {
  // Return the VD digital voltage in V:
  return (5.23);
}


void hal::setTBia(double /*IA*/) {
}

void hal::setTBva(double /*VA*/){
}

void hal::setTBid(double /*ID*/) {
}

void hal::setTBvd(double /*VD*/) {
}


bool hal::rocSetDACs(uint8_t /*rocId*/, std::map< uint8_t, uint8_t > dacPairs) {

  bool is_wbc = false;
  // Iterate over all DAC id/value pairs and set the DAC
  for(std::map< uint8_t,uint8_t >::iterator it = dacPairs.begin(); it != dacPairs.end(); ++it) {
    LOG(logDEBUGHAL) << "Set DAC" << static_cast<int>(it->first) << " to " << static_cast<int>(it->second);
    if(it->first == ROC_DAC_WBC) { is_wbc = true; }
  }

  // Make sure to issue a ROC Reset after WBC has been programmed:
  if(is_wbc) {
    LOG(logDEBUGHAL) << "WBC has been programmed - sending a ROC Reset command.";
    daqTriggerSingleSignal(TRG_SEND_RSR);
  }

  // Everything went all right:
  return true;
}

bool hal::rocSetDAC(uint8_t roci2c, uint8_t dacId, uint8_t dacValue) {

  LOG(logDEBUGHAL) << "ROC@I2C " << static_cast<size_t>(roci2c) 
		   << ": Set DAC" << static_cast<int>(dacId) << " to " << static_cast<int>(dacValue);

  // Make sure to issue a ROC Reset after the DAc WBC has been programmed:
  if(dacId == ROC_DAC_WBC) {
    LOG(logDEBUGHAL) << "WBC has been programmed - sending a ROC Reset command.";
    daqTriggerSingleSignal(TRG_SEND_RSR);
  }
  return true;
}

bool hal::tbmSetRegs(std::map< uint8_t, uint8_t > /*regPairs*/) {
  // Everything went all right:
  return true;
}

bool hal::tbmSetReg(uint8_t /*regId*/, uint8_t /*regValue*/) {
  return true;
}

void hal::RocSetMask(uint8_t /*rocid*/, bool /*mask*/, std::vector<pixelConfig> /*pixels*/) {
}

void hal::AllColumnsSetEnable(uint8_t /*rocid*/, bool /*enable*/) {
}

void hal::PixelSetCalibrate(uint8_t /*rocid*/, uint8_t /*column*/, uint8_t /*row*/, uint16_t /*flags*/) {
}

void hal::RocClearCalibrate(uint8_t /*rocid*/) {
}

void hal::SetupTrimValues(unsigned char, std::vector<pxar::pixelConfig, std::allocator<pxar::pixelConfig> >) {
}

void hal::SetupI2CValues(std::vector<unsigned char, std::allocator<unsigned char> >) {
}

// ---------------- TEST FUNCTIONS ----------------------

std::vector<Event*> hal::MultiRocAllPixelsCalibrate(std::vector<uint8_t> rocids, std::vector<int32_t> parameter) {

  uint32_t flags = static_cast<uint32_t>(parameter.at(0));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(1));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << nTriggers*ROC_NUMROWS*ROC_NUMCOLS << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t i = 0; i < ROC_NUMCOLS; i++) {
    for(size_t j = 0; j < ROC_NUMROWS; j++) {
      for(size_t k = 0; k < nTriggers; k++) {
	// New event:
	Event* evt = new Event();
	for(std::vector<uint8_t>::iterator roc = rocids.begin(); roc != rocids.end(); ++roc) {
	  // Fill the event with some pixels from the current ROC:
	  fillEvent(evt,*roc,i,j,flags);
	}
	// Count pixels:
	total_pixel += evt->pixels.size();
	data.push_back(evt);
	LOG(logDEBUGPIPES) << *evt;
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::MultiRocOnePixelCalibrate(std::vector<uint8_t> rocids, uint8_t column, uint8_t row, std::vector<int32_t> parameter) {

  uint32_t flags = static_cast<uint32_t>(parameter.at(0));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(1));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << nTriggers << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t k = 0; k < nTriggers; k++) {
    Event* evt = new Event();
    for(std::vector<uint8_t>::iterator roc = rocids.begin(); roc != rocids.end(); ++roc) {
      // Fill the event with some pixels from the current ROC:
      fillEvent(evt,*roc,column,row,flags);
    }
    // Count pixels:
    total_pixel += evt->pixels.size();
    data.push_back(evt);
    LOG(logDEBUGPIPES) << *evt;
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::SingleRocAllPixelsCalibrate(uint8_t rocid, std::vector<int32_t> parameter) {

  uint32_t flags = static_cast<uint32_t>(parameter.at(0));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(1));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << nTriggers*ROC_NUMROWS*ROC_NUMCOLS << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t i = 0; i < ROC_NUMCOLS; i++) {
    for(size_t j = 0; j < ROC_NUMROWS; j++) {
      for(size_t k = 0; k < nTriggers; k++) {
	// Create a new event:
	Event* evt = new Event();
	fillEvent(evt,rocid,i,j,flags);

	// Count pixels:
	total_pixel += evt->pixels.size();
	data.push_back(evt);
	LOG(logDEBUGPIPES) << *evt;
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::SingleRocOnePixelCalibrate(uint8_t rocid, uint8_t column, uint8_t row, std::vector<int32_t> parameter) {

  uint32_t flags = static_cast<uint32_t>(parameter.at(0));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(1));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << nTriggers << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t k = 0; k < nTriggers; k++) {
    Event* evt = new Event();
    fillEvent(evt,rocid,column,row,flags);

    // Count pixels:
    total_pixel += evt->pixels.size();
    data.push_back(evt);
    LOG(logDEBUGPIPES) << *evt;
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}


std::vector<Event*> hal::MultiRocAllPixelsDacScan(std::vector<uint8_t> rocids, std::vector<int32_t> parameter) {

  uint8_t dacmin = static_cast<uint8_t>(parameter.at(1));
  uint8_t dacmax = static_cast<uint8_t>(parameter.at(2));
  uint16_t flags = static_cast<uint16_t>(parameter.at(3));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(4));
  uint8_t dacstep = static_cast<uint8_t>(parameter.at(5));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dacmax-dacmin)/dacstep+1)*nTriggers*ROC_NUMROWS*ROC_NUMCOLS << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;
  uint8_t dachalf = static_cast<uint8_t>(dacmax-dacmin)/2;

  for(size_t i = 0; i < ROC_NUMCOLS; i++) {
    for(size_t j = 0; j < ROC_NUMROWS; j++) {
      for(size_t dac = 0; dac < static_cast<size_t>(dacmax-dacmin+1); dac += dacstep) {
	for(size_t k = 0; k < nTriggers; k++) {
	  // Create new event:
	  Event* evt = new Event();
	  for(std::vector<uint8_t>::iterator roc = rocids.begin(); roc != rocids.end(); ++roc) {
	    // Mimic some edge at 50% of the DAC range:
	    if((flags&FLAG_RISING_EDGE) && dac > dachalf) fillEvent(evt,*roc,i,j,flags);
	    else if(!(flags&FLAG_RISING_EDGE) && dac < dachalf) fillEvent(evt,*roc,i,j,flags);
	  }
	  // Count pixels:
	  total_pixel += evt->pixels.size();
	  data.push_back(evt);
	  LOG(logDEBUGPIPES) << *evt;
	}
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::MultiRocOnePixelDacScan(std::vector<uint8_t> rocids, uint8_t column, uint8_t row, std::vector<int32_t> parameter) {

  uint8_t dacmin = static_cast<uint8_t>(parameter.at(1));
  uint8_t dacmax = static_cast<uint8_t>(parameter.at(2));
  uint16_t flags = static_cast<uint16_t>(parameter.at(3));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(4));
  uint8_t dacstep = static_cast<uint8_t>(parameter.at(5));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dacmax-dacmin)/dacstep+1)*nTriggers << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;
  uint8_t dachalf = static_cast<uint8_t>(dacmax-dacmin)/2;

  for(size_t dac = 0; dac < static_cast<size_t>(dacmax-dacmin+1); dac += dacstep) {
    for(size_t k = 0; k < nTriggers; k++) {
      // Create new event:
      Event* evt = new Event();
      for(std::vector<uint8_t>::iterator roc = rocids.begin(); roc != rocids.end(); ++roc) {
	// Mimic some edge at 50% of the DAC range:
	if((flags&FLAG_RISING_EDGE) && dac > dachalf) fillEvent(evt,*roc,column,row,flags);
	else if(!(flags&FLAG_RISING_EDGE) && dac < dachalf) fillEvent(evt,*roc,column,row,flags);
      }
      // Count pixels:
      total_pixel += evt->pixels.size();
      data.push_back(evt);
      LOG(logDEBUGPIPES) << *evt;
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::SingleRocAllPixelsDacScan(uint8_t rocid, std::vector<int32_t> parameter) {

  uint8_t dacmin = static_cast<uint8_t>(parameter.at(1));
  uint8_t dacmax = static_cast<uint8_t>(parameter.at(2));
  uint16_t flags = static_cast<uint16_t>(parameter.at(3));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(4));
  uint8_t dacstep = static_cast<uint8_t>(parameter.at(5));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dacmax-dacmin)/dacstep+1)*nTriggers*ROC_NUMROWS*ROC_NUMCOLS << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;
  uint8_t dachalf = static_cast<uint8_t>(dacmax-dacmin)/2;

  for(size_t i = 0; i < ROC_NUMCOLS; i++) {
    for(size_t j = 0; j < ROC_NUMROWS; j++) {
      for(size_t dac = 0; dac < static_cast<size_t>(dacmax-dacmin+1); dac += dacstep) {
	for(size_t k = 0; k < nTriggers; k++) {
	  // Create new event:
	  Event* evt = new Event();
	  // Mimic some edge at 50% of the DAC range:
	  if((flags&FLAG_RISING_EDGE) && dac > dachalf) fillEvent(evt,rocid,i,j,flags);
	  else if(!(flags&FLAG_RISING_EDGE) && dac < dachalf) fillEvent(evt,rocid,i,j,flags);

	  // Count pixels:
	  total_pixel += evt->pixels.size();
	  data.push_back(evt);
	  LOG(logDEBUGPIPES) << *evt;
	}
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::SingleRocOnePixelDacScan(uint8_t rocid, uint8_t column, uint8_t row, std::vector<int32_t> parameter) {
  
  uint8_t dacmin = static_cast<uint8_t>(parameter.at(1));
  uint8_t dacmax = static_cast<uint8_t>(parameter.at(2));
  uint16_t flags = static_cast<uint16_t>(parameter.at(3));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(4));
  uint8_t dacstep = static_cast<uint8_t>(parameter.at(5));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dacmax-dacmin)/dacstep+1)*nTriggers << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;
  uint8_t dachalf = static_cast<uint8_t>(dacmax-dacmin)/2;

  for(size_t dac = 0; dac < static_cast<size_t>(dacmax-dacmin+1); dac += dacstep) {
    for(size_t k = 0; k < nTriggers; k++) {
      // Create a new event:
      Event* evt = new Event();
      // Mimic some edge at 50% of the DAC range:
      if((flags&FLAG_RISING_EDGE) && dac > dachalf) fillEvent(evt,rocid,column,row,flags);
      else if(!(flags&FLAG_RISING_EDGE) && dac < dachalf) fillEvent(evt,rocid,column,row,flags);

      // Count pixels:
      total_pixel += evt->pixels.size();
      data.push_back(evt);
      LOG(logDEBUGPIPES) << *evt;
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::MultiRocAllPixelsDacDacScan(std::vector<uint8_t> rocids, std::vector<int32_t> parameter) {

  uint8_t dac1min = static_cast<uint8_t>(parameter.at(1));
  uint8_t dac1max = static_cast<uint8_t>(parameter.at(2));
  uint8_t dac2min = static_cast<uint8_t>(parameter.at(4));
  uint8_t dac2max = static_cast<uint8_t>(parameter.at(5));
  uint16_t flags = static_cast<uint16_t>(parameter.at(6));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(7));
  uint8_t dac1step = static_cast<uint8_t>(parameter.at(8));
  uint8_t dac2step = static_cast<uint8_t>(parameter.at(9));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dac2max-dac2min)/dac2step+1)*static_cast<size_t>((dac1max-dac1min)/dac1step+1)*nTriggers*ROC_NUMROWS*ROC_NUMCOLS << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t i = 0; i < ROC_NUMCOLS; i++) {
    for(size_t j = 0; j < ROC_NUMROWS; j++) {
      for(size_t dac1 = 0; dac1 < static_cast<size_t>(dac1max-dac1min+1); dac1 += dac1step) {
	for(size_t dac2 = 0; dac2 < static_cast<size_t>(dac2max-dac2min+1); dac2 += dac2step) {
	  for(size_t k = 0; k < nTriggers; k++) {
	    // Create a new event:
	    Event* evt = new Event();
	    for(std::vector<uint8_t>::iterator roc = rocids.begin(); roc != rocids.end(); ++roc) {
	      // Mimic some working band of the two DACs:
	      if(isInTornadoRegion(dac1min, dac1max, dac1, dac2min, dac2max, dac2)) {
		fillEvent(evt,*roc,i,j,flags);
	      }
	    }
	    // Count pixels:
	    total_pixel += evt->pixels.size();
	    data.push_back(evt);
	    LOG(logDEBUGPIPES) << *evt;
	  }
	}
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::MultiRocOnePixelDacDacScan(std::vector<uint8_t> rocids, uint8_t column, uint8_t row, std::vector<int32_t> parameter) {

  uint8_t dac1min = static_cast<uint8_t>(parameter.at(1));
  uint8_t dac1max = static_cast<uint8_t>(parameter.at(2));
  uint8_t dac2min = static_cast<uint8_t>(parameter.at(4));
  uint8_t dac2max = static_cast<uint8_t>(parameter.at(5));
  uint16_t flags = static_cast<uint16_t>(parameter.at(6));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(7));
  uint8_t dac1step = static_cast<uint8_t>(parameter.at(8));
  uint8_t dac2step = static_cast<uint8_t>(parameter.at(9));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dac2max-dac2min)/dac2step+1)*static_cast<size_t>((dac1max-dac1min)/dac1step+1)*nTriggers << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t dac1 = 0; dac1 < static_cast<size_t>(dac1max-dac1min+1); dac1 += dac1step) {
    for(size_t dac2 = 0; dac2 < static_cast<size_t>(dac2max-dac2min+1); dac2 += dac2step) {
      for(size_t k = 0; k < nTriggers; k++) {
	// Create a new event:
	Event* evt = new Event();
	for(std::vector<uint8_t>::iterator roc = rocids.begin(); roc != rocids.end(); ++roc) {
	  // Mimic some working band of the two DACs:
	  if(isInTornadoRegion(dac1min, dac1max, dac1, dac2min, dac2max, dac2)) {
	    fillEvent(evt,*roc,column,row,flags);
	  }
	}
	// Count pixels:
	total_pixel += evt->pixels.size();
	data.push_back(evt);
	LOG(logDEBUGPIPES) << *evt;
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::SingleRocAllPixelsDacDacScan(uint8_t rocid, std::vector<int32_t> parameter) {

  uint8_t dac1min = static_cast<uint8_t>(parameter.at(1));
  uint8_t dac1max = static_cast<uint8_t>(parameter.at(2));
  uint8_t dac2min = static_cast<uint8_t>(parameter.at(4));
  uint8_t dac2max = static_cast<uint8_t>(parameter.at(5));
  uint16_t flags = static_cast<uint16_t>(parameter.at(6));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(7));
  uint8_t dac1step = static_cast<uint8_t>(parameter.at(8));
  uint8_t dac2step = static_cast<uint8_t>(parameter.at(9));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dac2max-dac2min)/dac2step+1)*static_cast<size_t>((dac1max-dac1min)/dac1step+1)*nTriggers*ROC_NUMROWS*ROC_NUMCOLS << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t i = 0; i < ROC_NUMCOLS; i++) {
    for(size_t j = 0; j < ROC_NUMROWS; j++) {
      for(size_t dac1 = 0; dac1 < static_cast<size_t>(dac1max-dac1min+1); dac1 += dac1step) {
	for(size_t dac2 = 0; dac2 < static_cast<size_t>(dac2max-dac2min+1); dac2 += dac2step) {
	  for(size_t k = 0; k < nTriggers; k++) {
	    // Create a new event:
	    Event* evt = new Event();
	    // Mimic some working band of the two DACs:
	    if(isInTornadoRegion(dac1min, dac1max, dac1, dac2min, dac2max, dac2)) {
	      fillEvent(evt,rocid,i,j,flags);
	    }
	    // Count pixels:
	    total_pixel += evt->pixels.size();
	    data.push_back(evt);
	    LOG(logDEBUGPIPES) << *evt;
	  }
	}
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

std::vector<Event*> hal::SingleRocOnePixelDacDacScan(uint8_t rocid, uint8_t column, uint8_t row, std::vector<int32_t> parameter) {

  uint8_t dac1min = static_cast<uint8_t>(parameter.at(1));
  uint8_t dac1max = static_cast<uint8_t>(parameter.at(2));
  uint8_t dac2min = static_cast<uint8_t>(parameter.at(4));
  uint8_t dac2max = static_cast<uint8_t>(parameter.at(5));
  uint16_t flags = static_cast<uint16_t>(parameter.at(6));
  uint16_t nTriggers = static_cast<uint16_t>(parameter.at(7));
  uint8_t dac1step = static_cast<uint8_t>(parameter.at(8));
  uint8_t dac2step = static_cast<uint8_t>(parameter.at(9));

  LOG(logDEBUGHAL) << "Flags: " << listFlags(flags);
  LOG(logDEBUGHAL) << "Expecting " << static_cast<size_t>((dac2max-dac2min)/dac2step+1)*static_cast<size_t>((dac1max-dac1min)/dac1step+1)*nTriggers << " events.";
  std::vector<Event*> data;
  size_t total_pixel = 0;

  for(size_t dac1 = 0; dac1 < static_cast<size_t>(dac1max-dac1min+1); dac1 += dac1step) {
    for(size_t dac2 = 0; dac2 < static_cast<size_t>(dac2max-dac2min+1); dac2 += dac2step) {
      for(size_t k = 0; k < nTriggers; k++) {
	// Create a new event:
	Event* evt = new Event();
	// Mimic some working band of the two DACs:
	if(isInTornadoRegion(dac1min, dac1max, dac1, dac2min, dac2max, dac2)) {
	  fillEvent(evt,rocid,column,row,flags);
	}
	// Count pixels:
	total_pixel += evt->pixels.size();
	data.push_back(evt);
	LOG(logDEBUGPIPES) << *evt;
      }
    }
  }

  LOG(logDEBUGHAL) << "Readout size: " << data.size() << " Events (" << total_pixel << " pixels).";
  return data;
}

// Testboard power switches:

void hal::HVon() {
  // Wait a little and let the HV relais do its job:
  mDelay(400);
}

void hal::HVoff() {
}
 
void hal::Pon() {
  // Wait a little and let the power switch do its job:
  mDelay(300);
}

void hal::Poff() {
}



// Testboard probe channel selectors:

void hal::SignalProbeD1(uint8_t /*signal*/) {
}

void hal::SignalProbeD2(uint8_t /*signal*/) {
}

void hal::SignalProbeA1(uint8_t /*signal*/) {
}

void hal::SignalProbeA2(uint8_t /*signal*/) {
}

void hal::SignalProbeADC(uint8_t /*signal*/, uint8_t /*gain*/) {
}

void hal::SetClockSource(uint8_t /*src*/) {
}

bool hal::IsClockPresent() {
  return true;
}

void hal::SetClockStretch(uint8_t /*src*/, uint16_t /*delay*/, uint16_t /*width*/) {
}

void hal::daqStart(uint8_t deser160phase, uint32_t buffersize) {

  LOG(logDEBUGHAL) << "Starting new DAQ session.";

  // Length of a token chain (number of ROCs per data stream):
  uint8_t tokenChainLength = 1; // One ROC for DESER160 readout.
  if(m_tbmtype != TBM_NONE && m_tbmtype != TBM_EMU) {
    // Four ROCs per stream for dual-400MHz, eight ROCs for single-400MHz readout:
    tokenChainLength *= (m_tbmtype >= TBM_09 ? 4 : 8);
    // Split the total buffer size when having more than one channel
    buffersize /= (m_tbmtype >= TBM_09 ? 4 : 2);
  }
  LOG(logDEBUGHAL) << "Determined Token Chain Length: " << static_cast<int>(tokenChainLength) << " ROCs.";

  uint32_t allocated_buffer_ch0 = buffersize;
  LOG(logDEBUGHAL) << "Allocated buffer size, Channel 0: " << allocated_buffer_ch0;
  src0 = dtbSource(NULL,0,tokenChainLength,m_tbmtype,m_roctype,true);
  src0 >> splitter0;

  if(m_tbmtype != TBM_NONE && m_tbmtype != TBM_EMU) {
    LOG(logDEBUGHAL) << "Enabling Deserializer400 for data acquisition.";

    uint32_t allocated_buffer_ch1 = buffersize;
    LOG(logDEBUGHAL) << "Allocated buffer size, Channel 1: " << allocated_buffer_ch1;
    src1 = dtbSource(NULL,1,tokenChainLength,m_tbmtype,m_roctype,true);
    src1 >> splitter1;

    // If we have an old TBM version set up the DESER400 to read old data format:
    // "old" is everything before TBM08B (so: TBM08, TBM08A)
    if(m_tbmtype < TBM_08B) { 
      LOG(logDEBUGHAL) << "Pre-series TBM with outdated trailer format. Configuring DESER400 accordingly.";
    }

    // For Dual-link TBMs (2x400MHz) we need even more DAQ channels:
    if(m_tbmtype >= TBM_09) {
      LOG(logDEBUGHAL) << "Dual-link TBM detected, enabling more DAQ channels.";

      uint32_t allocated_buffer_ch2 = buffersize;
      LOG(logDEBUGHAL) << "Allocated buffer size, Channel 2: " << allocated_buffer_ch2;
      src2 = dtbSource(NULL,2,tokenChainLength,m_tbmtype,m_roctype,true);
      src2 >> splitter2;

      uint32_t allocated_buffer_ch3 = buffersize;
      LOG(logDEBUGHAL) << "Allocated buffer size, Channel 3: " << allocated_buffer_ch3;
      src3 = dtbSource(NULL,3,tokenChainLength,m_tbmtype,m_roctype,true);
      src3 >> splitter3;
    }
  }
  else {
    LOG(logDEBUGHAL) << "Enabling Deserializer160 for data acquisition."
		     << " Phase: " << static_cast<int>(deser160phase);
  }
}

Event* hal::daqEvent() {

  Event* current_Event = new Event();

  return current_Event;
}

std::vector<Event*> hal::daqAllEvents() {

  std::vector<Event*> evt;
  return evt;
}

rawEvent* hal::daqRawEvent() {

  rawEvent* current_Event = new rawEvent();
  return current_Event;
}

std::vector<rawEvent*> hal::daqAllRawEvents() {

  std::vector<rawEvent*> raw;
  return raw;
}

std::vector<uint16_t> hal::daqBuffer() {
  
  std::vector<uint16_t> raw;
  return raw;
}

void hal::daqTriggerSource(uint16_t /*source*/) {}

void hal::daqTriggerSingleSignal(uint8_t /*signal*/) {

  // Attach the single signal direct source for triggers
  // in addition to the currently active source:
  LOG(logDEBUGHAL) << std::hex << TRG_SEL_SINGLE_DIR << " " << _currentTrgSrc
		   << " - " << (TRG_SEL_SINGLE_DIR | _currentTrgSrc) << std::dec;

}

void hal::daqTrigger(uint32_t /*nTrig*/, uint16_t /*period*/) {}

void hal::daqTriggerLoop(uint16_t /*period*/) {}

void hal::daqTriggerLoopHalt() {}

uint32_t hal::daqBufferStatus() { return 0; }

void hal::daqStop() {}

void hal::daqClear() {}

std::vector<uint16_t> hal::daqADC(uint8_t, uint8_t, uint16_t, uint8_t, uint8_t, uint8_t){
  return vector<uint16_t>();
}

void hal::SigSetMode(unsigned char, unsigned char) {}

std::vector<std::vector<uint16_t> > hal::daqReadback() {
  return std::vector<std::vector<uint16_t> >();
}

void hal::SigSetLVDS() {}

void hal::setHubId(unsigned char) {}

void hal::SigSetLCDS() {}

statistics hal::daqStatistics() {
  return statistics();
}

uint16_t hal::GetADC(uint8_t rpc_par1) {
  return rpc_par1;
}
