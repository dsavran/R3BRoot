#include "R3BBunchedFiberCal2Hit.h"
#include <cassert>
#include <TClonesArray.h>
#include <FairLogger.h>
#include "R3BTCalEngine.h"
#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"

R3BBunchedFiberCal2Hit::ToT::ToT(R3BBunchedFiberCalData const *a_lead,
    R3BBunchedFiberCalData const *a_trail, Double_t a_tot)
  : lead(a_lead)
  , trail(a_trail)
  , tot(a_tot)
{
}

R3BBunchedFiberCal2Hit::R3BBunchedFiberCal2Hit(const char *a_name, Int_t
    a_verbose, Direction a_direction, UInt_t a_sub_num, UInt_t a_mapmt_per_sub,
    UInt_t a_spmt_per_sub)
  : FairTask(TString("R3B") + a_name + "Cal2Hit", a_verbose)
  , fName(a_name)
  , fDirection(a_direction)
  , fSubNum(a_sub_num)
  , fCalItems(nullptr)
  , fHitItems(new TClonesArray("R3BBunchedFiberHitData"))
  , fNofHitItems(0)
{
  fChPerSub[0] = a_mapmt_per_sub;
  fChPerSub[1] = a_spmt_per_sub;
}

R3BBunchedFiberCal2Hit::~R3BBunchedFiberCal2Hit()
{
  delete fHitItems;
}

InitStatus R3BBunchedFiberCal2Hit::Init()
{
  auto mgr = FairRootManager::Instance();
  if (!mgr) {
    LOG(ERROR) << "FairRootManager not found." << FairLogger::endl;
    return kERROR;
  }
  auto name = fName + "Cal";
  fCalItems = (TClonesArray *)mgr->GetObject(name);
  if (!fCalItems) {
    LOG(ERROR) << "Branch " << name << " not found." << FairLogger::endl;
    return kERROR;
  }
  mgr->Register(fName + "Hit", "Land", fHitItems, kTRUE);
  // Resize per-channel info arrays.
  for (auto side_i = 0; side_i < 2; ++side_i) {
    fChannelArray[side_i].resize(fSubNum * fChPerSub[side_i]);
  }
  return kSUCCESS;
}

InitStatus R3BBunchedFiberCal2Hit::ReInit()
{
  return kSUCCESS;
}

void R3BBunchedFiberCal2Hit::Exec(Option_t *option)
{
  for (auto side_i = 0; side_i < 2; ++side_i) {
    // Clear local helper containers.
    auto &array = fChannelArray[side_i];
    for (auto it = array.begin(); array.end() != it; ++it) {
      it->lead_list.clear();
      it->tot_list.clear();
    }
  }

  auto const c_period = 4096e3 / CLOCK_TDC_MHZ;
  size_t cal_num = fCalItems->GetEntriesFast();

  // Find multi-hit ToT for every channel.
  // The easiest safe way to survive ugly cases is to record all
  // leading edges per channel, and then pair up with whatever
  // trailing we have.
  // Not super efficient, but shouldn't crash if the data is not
  // perfect.
  unsigned n_lead = 0;
  unsigned n_trail = 0;
  for (size_t j = 0; j < cal_num; ++j) {
    auto cur_cal = (R3BBunchedFiberCalData const *)fCalItems->At(j);
    if (cur_cal->IsLeading()) {
	  ++n_lead;
      auto side_i = cur_cal->IsMAPMT() ? 0 : 1;
	  auto ch_i = cur_cal->GetChannel();
	  auto &channel = fChannelArray[side_i].at(ch_i - 1);
	  channel.lead_list.push_back(cur_cal);
    } else {
	  ++n_trail;
	}
  }
  if (n_lead != n_trail) {
	  return;
  }
  for (size_t j = 0; j < cal_num; ++j) {
    auto cur_cal = (R3BBunchedFiberCalData const *)fCalItems->At(j);
    if (cur_cal->IsTrailing()) {
	  auto side_i = cur_cal->IsMAPMT() ? 0 : 1;
	  auto ch_i = cur_cal->GetChannel();
      auto &channel = fChannelArray[side_i].at(ch_i - 1);
      if (channel.lead_list.empty()) {
		break;
	  }
      auto lead = channel.lead_list.front();
      channel.lead_list.pop_front();
      auto tot = fmod(cur_cal->GetTime_ns() - lead->GetTime_ns() + c_period, c_period);
      channel.tot_list.push_back(ToT(lead, cur_cal, tot));
    }
  }

  // Make every permutation to create fibers.
  // TODO: Some combinations can probably be excluded, don't know yet how.
  auto const &mapmt_array = fChannelArray[0];
  auto const &spmt_array = fChannelArray[1];
  for (auto it_mapmt = mapmt_array.begin(); mapmt_array.end() != it_mapmt;
      ++it_mapmt) {
    auto const &mapmt = *it_mapmt;
    for (auto it_mapmt_tot = mapmt.tot_list.begin(); mapmt.tot_list.end() !=
        it_mapmt_tot; ++it_mapmt_tot) {
      auto const &mapmt_tot = *it_mapmt_tot;
      for (auto it_spmt = spmt_array.begin(); spmt_array.end() != it_spmt; ++it_spmt) {
        auto const &spmt = *it_spmt;
        for (auto it_spmt_tot = spmt.tot_list.begin(); spmt.tot_list.end() !=
            it_spmt_tot; ++it_spmt_tot) {
          auto const &spmt_tot = *it_spmt_tot;

          /*
           * How to get a fiber ID for a fiber detector defined as:
           *  SubNum = 2
           *  MAPMT = 256
           *  SPMT = 2
           * This means we'll receive 512 MAPMT channels as 1..512, and 4 SPMT
           * channels, but the first half (sub-detector) is completely
           * decoupled from the second half. The sequence of all fibers in
           * order is then, as (MAPMT,SPMT)-pairs:
           *  (1,1), (1,2), (2,1), ... (256,1), (256,2),
           *  (257,3), (257,4), (258,3), ... (512,3), (512,4)
           */
          auto fiber_id = (mapmt_tot.lead->GetChannel() - 1) * fChPerSub[1] +
            (spmt_tot.lead->GetChannel() % fChPerSub[1]);

          // TODO: Use it_sub->direction to find real life coordinates.

          // Fix fiber installation mistakes.
          fiber_id = FixMistake(fiber_id);

          new ((*fHitItems)[fNofHitItems++])
            R3BBunchedFiberHitData(fiber_id,
                mapmt_tot.lead->GetTime_ns(),
                spmt_tot.lead->GetTime_ns(),
                mapmt_tot.tot,
                spmt_tot.tot);
        }
      }
    }
  }
}

void R3BBunchedFiberCal2Hit::FinishEvent()
{
  fHitItems->Clear();
  fNofHitItems = 0;
}

void R3BBunchedFiberCal2Hit::FinishTask()
{
}

ClassImp(R3BBunchedFiberCal2Hit)
