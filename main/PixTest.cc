#include <iostream>
#include "PixTest.hh"

using namespace std;

ClassImp(PixTest)

// ----------------------------------------------------------------------
PixTest::PixTest(TBInterface *tb, string name, PixTestParameters *tp) {
  init(tb, name, tp); 
  cout << "PixTest ctor(TBInterface *, string)" << endl;
}


// ----------------------------------------------------------------------
PixTest::PixTest() {
  cout << "PixTestBase ctor()" << endl;

}


// ----------------------------------------------------------------------
void PixTest::init(TBInterface *tb, string name, PixTestParameters *tp) {
  cout << "PixTest::init()" << endl;
  fTB = tb; 
  fName = name;
  fParameters = tp->getTestParameters(name); 
  NCOL = 52; 
  NROW = 80;

  for (map<string,string>::iterator imap = fParameters.begin(); imap != fParameters.end(); ++imap) {
    setParameter(imap->first, imap->second); 
  }
}


// // ----------------------------------------------------------------------
// bool PixTest::getParameter(string parName) {
//   map<string, string> testPars = fTestParameters->getTestParameters(fName); 
//   return false; 
// }

// ----------------------------------------------------------------------
bool PixTest::setParameter(string parName, string value) {
  cout << " PixTest::setParameter wrong function" << endl;
  return false;
}


// ----------------------------------------------------------------------
bool PixTest::getParameter(std::string parName, int &ival) {
  bool found(false);
  for (map<string,string>::iterator imap = fParameters.begin(); imap != fParameters.end(); ++imap) {
    if (0 == imap->first.compare(parName)) {
      found = true; 
      break;
    }
  }
  if (found) {
    ival = atoi(fParameters[parName].c_str()); 
  }
  return found; 
}


  // ----------------------------------------------------------------------
bool PixTest::getParameter(std::string parName, float &fval) {
  bool found(false);
  for (map<string,string>::iterator imap = fParameters.begin(); imap != fParameters.end(); ++imap) {
    if (0 == imap->first.compare(parName)) {
      found = true; 
      break;
    }
  }
  if (found) {
    fval = atof(fParameters[parName].c_str()); 
  }
  return found; 
}


// ----------------------------------------------------------------------
void PixTest::dumpParameters() {
  cout << "Parameters for test" << getName() << endl;
  for (map<string,string>::iterator imap = fParameters.begin(); imap != fParameters.end(); ++imap) {
    cout << imap->first << ": " << imap->second << endl;
  }
}


// ----------------------------------------------------------------------
PixTest::~PixTest() {
  cout << "PixTestBase dtor()" << endl;
}

// ----------------------------------------------------------------------
void PixTest::testDone() {
  cout << "PixTest::testDone()" << endl;
  Emit("testDone()"); 
}

// ----------------------------------------------------------------------
void PixTest::update() {
  //  cout << "PixTest::update()" << endl;
  Emit("update()"); 
}


// ----------------------------------------------------------------------
void PixTest::doTest() {
  cout << "PixTest::doTest()" << endl;
}


// ----------------------------------------------------------------------
void PixTest::doAnalysis() {
  cout << "PixTest::doAnalysis()" << endl;

}


// ----------------------------------------------------------------------
TH1* PixTest::nextHist() {
  cout << "PixTest::nextHist()" << endl;
//   if (fDisplayedHist) {
//     TH1 *h(0); 
//     for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
//       if (fDisplayedHist == (*il)) {
// 	// missing "back()"
// 	h = (*il); 
// 	fDisplayedHist = h; 
// 	return h; 
//       }
//     }
//   }
  return 0; 
}

// ----------------------------------------------------------------------
TH1* PixTest::previousHist() {
  cout << "PixTest::previousHist()" << endl;
//   if (fDisplayedHist) {
//     TH1 *h(0); 
//     for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
//       if (fDisplayedHist == (*il)) {
// 	// missing "back()"
// 	h = *il; 
// 	fDisplayedHist = h; 
// 	return h; 
//       }
//     }
//   }
  return 0; 
}

// ----------------------------------------------------------------------
void PixTest::setTitles(TH1 *h, const char *sx, const char *sy, float size, 
			float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy); 
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

// ----------------------------------------------------------------------
void PixTest::clearHist() {
  for (list<TH1*>::iterator il = fHistList.begin(); il != fHistList.end(); ++il) {
    (*il)->Reset();
  }
}