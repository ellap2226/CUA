#ifndef analysis_TestAnalysis_h
#define analysis_TestAnalysis_h

#include "SampleAnalyzer/Process/Analyzer/AnalyzerBase.h"
#include <TH1F.h>

namespace MA5
{
  class TestAnalysis : public AnalyzerBase
  {
    INIT_ANALYSIS(TestAnalysis,"TestAnalysis")

    public:
    virtual bool Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters);
    virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
    virtual bool Execute(SampleFormat& sample, const EventFormat& event);

  private:
    TH1F* ptaamtotal;
    TH1F* ptbbmtotal;
    TH1F* ptamaal;
    TH1F* ptamaasl;
    TH1F* ptbmbbl;
    TH1F* ptbmbbsl;
    TH1F* minAngularDist;
    TH1F* cosThetaY;
    TH1F* cosThetaH;
    TH1F* cosThetaPho1;
    TH1F* cosThetaPho2;
    TH1F* cosThetaB;
    TH1F* cosThetaAB;
    /*
    ptaamtotal->Draw();
    ptbbmtotal->Draw();
    ptamaal->Draw();
    ptamaasl->Draw();
    ptbmbbl->Draw();
    ptbmbbsl->Draw();
    minAngularDist->Draw();
    cosThetaY->Draw();
    cosThetaH->Draw();
    cosThetaPho1->Draw();
    cosThetaPho2->Draw();
    cosThetaB->Draw();
    cosThetaAB->Draw();
    */
    
  };
}

#endif
