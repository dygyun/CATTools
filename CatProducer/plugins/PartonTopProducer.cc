#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace edm;
using namespace reco;

class PartonTopProducer : public edm::EDProducer
{
public:
  PartonTopProducer(const edm::ParameterSet& pset);
  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

  enum TTbarMode { CH_NONE = 0, CH_FULLHADRON = 1, CH_SEMILEPTON, CH_FULLLEPTON };
  enum DecayMode { CH_HADRON = 1, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

private:
  const reco::Candidate* getLast(const reco::Candidate* p) const;
  reco::GenParticleRef buildGenParticle(const reco::Candidate* p, reco::GenParticleRefProd& refHandle,
                                        std::auto_ptr<reco::GenParticleCollection>& outColl) const;

  typedef reco::Particle::LorentzVector LorentzVector;

private:
  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken_;
};

PartonTopProducer::PartonTopProducer(const edm::ParameterSet& pset)
{
  genParticleToken_ = consumes<edm::View<reco::Candidate> >(pset.getParameter<edm::InputTag>("genParticles"));

  produces<reco::GenParticleCollection>();
  produces<int>("channel");
  produces<std::vector<int> >("modes");
}

void PartonTopProducer::produce(edm::Event& event, const edm::EventSetup& eventSetup)
{
  edm::Handle<edm::View<reco::Candidate> > genParticleHandle;
  event.getByToken(genParticleToken_, genParticleHandle);

  std::auto_ptr<reco::GenParticleCollection> partons(new reco::GenParticleCollection);
  auto partonRefHandle = event.getRefBeforePut<reco::GenParticleCollection>();

  std::auto_ptr<int> channel(new int(CH_NONE));
  std::auto_ptr<std::vector<int> > modes(new std::vector<int>());

  // Collect top quarks and unstable B-hadrons
  std::vector<const reco::Candidate*> tQuarks;
  for ( size_t i=0, n=genParticleHandle->size(); i<n; ++i )
  {
    const reco::Candidate& p = genParticleHandle->at(i);
    const int status = p.status();
    if ( status == 1 ) continue;

    // Collect parton level objects. Ignore gluons and photons
    const int absPdgId = abs(p.pdgId());
    if ( absPdgId == 6 )
    {
      // top quark : select one 'after radiations'
      bool toKeep = true;
      if ( p.numberOfDaughters() == 0 ) toKeep = false;
      for ( size_t j=0, m=p.numberOfDaughters(); j<m; ++j )
      {
        const int dauId = p.daughter(j)->pdgId();
        if ( dauId == p.pdgId() ) { toKeep = false; break; }
      }
      if ( toKeep ) tQuarks.push_back(&p);
    }
  }
  // Build top quark decay tree in parton level
  // Also determine decay mode from parton level information
  size_t nElectron = 0, nMuon = 0, nTau = 0, nTauToLepton = 0;
  for ( int i=0, n=tQuarks.size(); i<n; ++i )
  {
    const reco::Candidate* t = tQuarks.at(i);
    const reco::Candidate* tLast = getLast(t);
    reco::GenParticleRef tRef = buildGenParticle(tLast, partonRefHandle, partons);
  }

  for ( int i=0, n=tQuarks.size(); i<n; ++i )
  {
    const reco::Candidate* tLast = getLast(tQuarks.at(i));
    reco::GenParticleRef tRef(partonRefHandle, i);

    const reco::Candidate* w = 0;
    const reco::Candidate* b = 0;
    for ( int j=0, m=tLast->numberOfDaughters(); j<m; ++j )
    {
      const reco::Candidate* dau = tLast->daughter(j);
      const unsigned int dauAbsId = abs(dau->pdgId());
      if ( dauAbsId == 24 and !w ) w = dau;
      else if ( dauAbsId < 6 and !b ) b = dau;
    }
    if ( !w or !b ) continue;
    reco::GenParticleRef wRef = buildGenParticle(w, partonRefHandle, partons);
    reco::GenParticleRef bRef = buildGenParticle(b, partonRefHandle, partons);
    partons->at(wRef.key()).addMother(tRef);
    partons->at(bRef.key()).addMother(tRef);
    partons->at(tRef.key()).addDaughter(wRef);
    partons->at(tRef.key()).addDaughter(bRef);

    // W decay products
    const reco::Candidate* wLast = getLast(w);
    const reco::Candidate* wDau1 = 0;
    const reco::Candidate* wDau2 = 0;
    for ( int j=0, m=wLast->numberOfDaughters(); j<m; ++j )
    {
      const reco::Candidate* dau = wLast->daughter(j);
      const unsigned int dauAbsId = abs(dau->pdgId());
      if ( dauAbsId > 16 ) continue; // Consider quarks and leptons only for W decays

      if ( !wDau1 ) wDau1 = dau;
      else if ( !wDau2 ) wDau2 = dau;
      else
      {
        cout << "--------------------------------------" << endl;
        cout << "WDECAY with more than 2 body!!! " << wLast->numberOfDaughters() << endl;
        cout << " dau1 = " << wDau1->pdgId() << " dau2 = " << wDau2->pdgId() << " skipped = " << dau->pdgId() << endl;
        cout << "--------------------------------------" << endl;
      }
    }
    if ( !wDau1 or !wDau2 ) continue;
    if ( abs(wDau1->pdgId()) > abs(wDau2->pdgId()) ) swap(wDau1, wDau2);
    reco::GenParticleRef wDauRef1 = buildGenParticle(wDau1, partonRefHandle, partons);
    reco::GenParticleRef wDauRef2 = buildGenParticle(wDau2, partonRefHandle, partons);
    partons->at(wDauRef1.key()).addMother(wRef);
    partons->at(wDauRef2.key()).addMother(wRef);
    partons->at(wRef.key()).addDaughter(wDauRef1);
    partons->at(wRef.key()).addDaughter(wDauRef2);

    // Special care for tau->lepton decays
    // Note : we do not keep neutrinos from tau decays (tau->W, nu_tau, W->l, nu_l)
    // Note : Up to 6 neutrinos from top decays if both W decays to taus and all taus go into leptonic decay chain
    const reco::Candidate* lepFromTau = 0;
    if ( abs(wDau1->pdgId()) == 15 )
    {
      const reco::Candidate* tauLast = getLast(wDau1);
      for ( int j=0, m=tauLast->numberOfDaughters(); j<m; ++j )
      {
        const reco::Candidate* dau = tauLast->daughter(j);
        const unsigned int dauAbsId = abs(dau->pdgId());
        if ( dauAbsId == 11 or dauAbsId == 13 )
        {
          if ( !lepFromTau ) lepFromTau = dau;
          else
          {
            cout << "--------------------------------------" << endl;
            cout << "TAU decay with more than 1 leptons!!!, nDau=" << tauLast->numberOfDaughters() << endl;
            cout << " dau = " << lepFromTau->pdgId() << " skipped = " << dau->pdgId() << endl;
            cout << "--------------------------------------" << endl;
          }
        }
      }
      if ( lepFromTau )
      {
        reco::GenParticleRef lepRef = buildGenParticle(lepFromTau, partonRefHandle, partons);
        partons->at(lepRef.key()).addMother(wDauRef1);
        partons->at(wDauRef1.key()).addDaughter(lepRef);
      }
    }
    int mode = 0;
    switch ( abs(wDau1->pdgId()) )
    {
      case 11: ++nElectron; mode = CH_ELECTRON; break;
      case 13: ++nMuon; mode = CH_MUON; break;
      case 15:
        ++nTau; mode = CH_TAU_HADRON;
        if ( lepFromTau )
        {
          ++nTauToLepton;
          if ( abs(lepFromTau->pdgId()) == 13 ) mode += 1;
          else mode += 2;
        }
        break;
    }
    modes->push_back(mode);
  }

  if ( modes->size() == 2 )
  {
    const int nLepton = nElectron + nMuon;
    if      ( nLepton == 0 ) *channel = CH_FULLHADRON;
    else if ( nLepton == 1 ) *channel = CH_SEMILEPTON;
    else if ( nLepton == 2 ) *channel = CH_FULLLEPTON;
  }

  event.put(partons);
  event.put(channel, "channel");
  event.put(modes, "modes");
}

const reco::Candidate* PartonTopProducer::getLast(const reco::Candidate* p) const
{
  int nDecay = 0;
  std::vector<const reco::Candidate*> sameCopies;
  for ( size_t i=0, n=p->numberOfDaughters(); i<n; ++i )
  {
    const reco::Candidate* dau = p->daughter(i);
    const int dauId = dau->pdgId();
    if ( dauId == 22 or dauId == 21 ) continue;
    if ( p->pdgId() == dau->pdgId() ) sameCopies.push_back(dau);
    else ++nDecay;
  }
  if ( nDecay == 0 )
  {
    for ( const auto dau : sameCopies )
    {
      if ( p->pdgId() == dau->pdgId() ) return getLast(dau);
    }
  }
  return p;
}

reco::GenParticleRef PartonTopProducer::buildGenParticle(const reco::Candidate* p, reco::GenParticleRefProd& refHandle,
                                                               std::auto_ptr<reco::GenParticleCollection>& outColl) const
{
  reco::GenParticle pOut(*dynamic_cast<const reco::GenParticle*>(p));
  pOut.clearMothers();
  pOut.clearDaughters();
  pOut.resetMothers(refHandle.id());
  pOut.resetDaughters(refHandle.id());

  outColl->push_back(pOut);

  return reco::GenParticleRef(refHandle, outColl->size()-1);
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PartonTopProducer);

