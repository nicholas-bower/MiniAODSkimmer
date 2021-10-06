// -*- C++ -*-
//
// Package:    MiniAODCleanerTest/MiniAODCleaner
// Class:      MediumLowPtElectronFilter
// 
/**\class MediumLowPtElectronFilter MediumLowPtElectronFilter.cc MiniAODCleanerTest/MiniAODCleaner/plugins/MediumLowPtElectronFilter.cc

 Description: [one line class summary]
Creates a collection of electrons passing the Loose Electron ID.returns true if at least one electron passes.
Recomissioned for use in MiniAOD cleaning.

Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Tue, 01 Jun 2021 
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "PhysicsTools/SelectorUtils/interface/CutApplicatorBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"


#include "TMath.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/RefToBaseVector.h"

using namespace edm;
using namespace std;
//
// class declaration
//

class MediumLowPtElectronFilter : public edm::stream::EDFilter<> {
public:
      explicit MediumLowPtElectronFilter(const edm::ParameterSet&);
  ~MediumLowPtElectronFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  
  
private:
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;
  edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MediumLowPtElectronFilter::MediumLowPtElectronFilter(const edm::ParameterSet& iConfig):
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("LowPtElectrons")))

  //trk_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("Tracks")))
{
  //now do what ever initialization is needed
  produces<pat::ElectronCollection>( "MiniMediumLowPtElectron" );
  produces<pat::ElectronRefVector>("MediumLowPtElectronRef");
}


MediumLowPtElectronFilter::~MediumLowPtElectronFilter()
{
  
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  
}

// ------------ method called on each new Event  ------------
bool MediumLowPtElectronFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace reco;

 
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_,electrons);
  unique_ptr<pat::ElectronCollection> passedelectrons(new pat::ElectronCollection);
  unique_ptr<pat::ElectronRefVector> passedelectronRef(new pat::ElectronRefVector);
  Handle<reco::VertexCollection> Vertex;


  for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end(); ++iele)
    { 
      pat::ElectronRef ERef(electrons,iele-electrons->begin()); 

	  if((iele->pt()>1)&&(iele->electronID("ID")>2))

	    {


	      passedelectrons->push_back(*iele);
	      passedelectronRef->push_back(ERef);
	      cout<< "EEPushback" <<endl;

	    }

	}


  iEvent.put(move(passedelectrons), "MiniMeduimElectron");
  iEvent.put(move(passedelectronRef),"MediumElectronRef");
  return true;
  //else return false;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
MediumLowPtElectronFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
MediumLowPtElectronFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
  void
  MediumLowPtElectronFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
MediumLowPtElectronFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MediumLowPtElectronFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MediumLowPtElectronFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MediumLowPtElectronFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MediumLowPtElectronFilter);
