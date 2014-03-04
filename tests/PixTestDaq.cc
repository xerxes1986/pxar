#include <stdlib.h>     /* atof, atoi */
#include <algorithm>    // std::find
#include <iostream>
#include "PixTestDaq.hh"
#include "log.h"

#include <TH2.h>

using namespace std;
using namespace pxar;

ClassImp(PixTestDaq)

// ----------------------------------------------------------------------
PixTestDaq::PixTestDaq(PixSetup *a, std::string name) : PixTest(a, name), fParNtrig(-1) {
  PixTest::init();
  init(); 
  LOG(logDEBUG) << "PixTestDaq ctor(PixSetup &a, string, TGTab *)";
}


//----------------------------------------------------------
PixTestDaq::PixTestDaq() : PixTest() {
  LOG(logDEBUG) << "PixTestDaq ctor()";
}

// ----------------------------------------------------------------------
bool PixTestDaq::setParameter(string parName, string sval) {
  bool found(false);
  string stripParName; 
  for (unsigned int i = 0; i < fParameters.size(); ++i) {
    if (fParameters[i].first == parName) {
      found = true; 
      
      LOG(logDEBUG) << "  ==> parName: " << parName;
      LOG(logDEBUG) << "  ==> sval:    " << sval;
      if (!parName.compare("Ntrig")) {
	fParNtrig = atoi(sval.c_str()); 
	setToolTips();
      }
      break;
    }
  }
  return found; 
}


// ----------------------------------------------------------------------
void PixTestDaq::init() {
  LOG(logDEBUG) << "PixTestDaq::init()";

  setToolTips();
  fDirectory = gFile->GetDirectory(fName.c_str()); 
  if (!fDirectory) {
    fDirectory = gFile->mkdir(fName.c_str()); 
  } 
  fDirectory->cd(); 

}

// ----------------------------------------------------------------------
void PixTestDaq::setToolTips() {
  fTestTip    = string("run DAQ")
    ;
  fSummaryTip = string("to be implemented")
    ;
}


// ----------------------------------------------------------------------
void PixTestDaq::bookHist(string name) {
  fDirectory->cd(); 
  LOG(logDEBUG) << "nothing done with " << name; 
}


//----------------------------------------------------------
PixTestDaq::~PixTestDaq() {
  LOG(logDEBUG) << "PixTestDaq dtor";
}


// ----------------------------------------------------------------------
void PixTestDaq::doTest() {

  LOG(logINFO) << "PixTestDaq::doTest() start with fParNtrig = " << fParNtrig;

  PixTest::update(); 
  fDirectory->cd();
  vector<TH2D*> hits;


  // All on!
  fApi->_dut->testAllPixels(false);
  //fApi->_dut->maskAllPixels(true);

  // Set some pixels up for getting calibrate signals:
  for (int i = 0; i < 3; ++i) {
    fApi->_dut->testPixel(i, 5, true);
    fApi->_dut->maskPixel(i, 5, false);
    fApi->_dut->testPixel(i, 6, true);
    fApi->_dut->maskPixel(i, 6, false);
  }

    
  // Start the DAQ:
  fApi->daqStart(fPixSetup->getConfigParameters()->getTbPgSettings());
  
  // Send the triggers:
  fApi->daqTrigger(fParNtrig);
   
  // Stop the DAQ:
  fApi->daqStop();

  // And read out the full buffer:
  // Either fetching undecoded raw data events or decoded events, check API documentation!
  vector<pxar::Event> daqdat = fApi->daqGetEventBuffer();

  cout << "Number of events read from board: " << daqdat.size() << endl;

  string name("Hits");
  vector<uint8_t> rocIds = fApi->_dut->getEnabledRocIDs();
  TH2D *h2(0);

  for (unsigned int iroc = 0; iroc < rocIds.size(); ++iroc){
    id2idx.insert(make_pair(rocIds[iroc], iroc));
    h2 = bookTH2D(Form("%s_C%d", name.c_str(), iroc), Form("%s_C%d", name.c_str(), rocIds[iroc]), 52, 0., 52., 80, 0., 80.);
    h2->SetMinimum(0.);
    h2->SetDirectory(fDirectory);
    setTitles(h2, "col", "row");
    hits.push_back(h2);
  }

 
  for(std::vector<pxar::Event>::iterator it = daqdat.begin(); it != daqdat.end(); ++it) {
    LOG(logDEBUG) << (*it) << std::endl;

   
    hits[id2idx[it->roc_id]]->SetBinContent(it->column,it->row,it->value);
 
  }
  // We should store the data, probably...


  for (unsigned int i = 0; i < hits.size(); ++i) {
    fHistOptions.insert(make_pair(hits[i], "colz"));
  }

  copy(hitss.begin(), hits.end(), back_inserter(fHistList));


  TH2D *h = (TH2D*)(*fHistList.begin());

  h->Draw(getHistOption(h).c_str());
  fDisplayedHist = find(fHistList.begin(), fHistList.end(), h);

  PixTest::update();
 

  LOG(logINFO) << "PixTestDaq::doTest() done";
}
