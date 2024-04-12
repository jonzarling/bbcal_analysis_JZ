#include <map>
using namespace std;

#include "DEventProcessor_BBCAL_tree.h"

#include <DANA/DApplication.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseData.h>
#include <PAIR_SPECTROMETER/DPSPair.h>
#include <PAIR_SPECTROMETER/DPSCPair.h>
#include <TRIGGER/DL1Trigger.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGMHit.h>


// Routine used to create our DEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new DEventProcessor_BBCAL_tree());
    
  }
} // "C"


// Returns 0,channel for north
// Returns 1,channel for south
pair<int,int> DEventProcessor_BBCAL_tree::getBBCALIndex(int crate, int slot, int channel) {
	// For EIC-BCAL runs
	// crate : 94, slot 7 for North and 8 for South and Channels from 0 to 15.
	if(crate!=94 || slot<7 || slot>8 || channel<0 || channel>15) return make_pair(-1,-1);
	return make_pair(slot-7,channel);
}


//------------------
// init
//------------------
jerror_t DEventProcessor_BBCAL_tree::init(void){
   	
	// Read in constants from gain file, if supplied
	string gainfile = "miniBCAL_constants.in";
	gPARMS->SetDefaultParameter("BBCAL:gainfile",gainfile);
	myfile = fopen(gainfile.c_str(), "r");
	
	// Check for flag to switch to cosmic configuration (default false)
	gPARMS->SetDefaultParameter("BBCAL:isCosmicSetup",isCosmicSetup);
	
	// Check for flag to switch to cosmic configuration (default false)
	gPARMS->SetDefaultParameter("BBCAL:isRawMode",isRawMode);
	
	// ROOT output filename
	string treename = "BBCAL_tree.root";
	gPARMS->SetDefaultParameter("BBCAL:FNAME",treename);

	//CREATE TTREE, TFILE
	dTreeInterface = DTreeInterface::Create_DTreeInterface("tr",treename); // "treename","filename"
	//TTREE BRANCHES
	DTreeBranchRegister locBranchRegister;
	locBranchRegister.Register_Single<int32_t>("run");
	locBranchRegister.Register_Single<ULong64_t>("event");
	
	if(!isCosmicSetup) {
		locBranchRegister.Register_Single<Int_t>("PS_NumPairs");
		locBranchRegister.Register_Single<Int_t>("PSC_NumPairs");
		locBranchRegister.Register_Single<Int_t>("PS_column_left");
		locBranchRegister.Register_Single<Int_t>("PS_column_right");
		locBranchRegister.Register_Single<Double_t>("PS_column_right_E");
		locBranchRegister.Register_Single<Double_t>("PS_column_left_E");
		
		// locBranchRegister.Register_Single<Int_t>("NTAGM_PASSING");
		// locBranchRegister.Register_FundamentalArray<Double_t>("Ediff_tagm","NTAGM_PASSING", 20); // Stringname, array_size stringname, init size
		// locBranchRegister.Register_FundamentalArray<Double_t>("tdiff_tagm","NTAGM_PASSING", 20); // Stringname, array_size stringname, init size
		// locBranchRegister.Register_Single<Int_t>("NTAGH_PASSING");
		// locBranchRegister.Register_FundamentalArray<Double_t>("Ediff_tagh","NTAGH_PASSING", 20); // Stringname, array_size stringname, init size
		// locBranchRegister.Register_FundamentalArray<Double_t>("tdiff_tagh","NTAGH_PASSING", 20); // Stringname, array_size stringname, init size
		
		// Integral reported in cosmic data (no PS trigger req.) would need calculating by hand. Save confusion by not saving a garbage number in cosmic data.
		for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("int_N"+to_string(i)); // Pulse integral  (pedestal subtracted, ADC units)
		for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("int_S"+to_string(i)); // Pulse integral  (pedestal subtracted, ADC units)
	}
	
	if(myfile!=NULL) {
		ParseGainFile(gainfile);
		locBranchRegister.Register_Single<Float_t>("energy_miniBCAL");
		locBranchRegister.Register_Single<Float_t>("energy_miniBCAL_2Sided");
		locBranchRegister.Register_Single<Float_t>("energy_miniBCAL_N");
		locBranchRegister.Register_Single<Float_t>("energy_miniBCAL_S");
		locBranchRegister.Register_Single<Float_t>("energy_miniBCAL_N_2Sided"); // Only uses channels if both N+S fire
		locBranchRegister.Register_Single<Float_t>("energy_miniBCAL_S_2Sided"); // Only uses channels if both N+S fire
		for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("energy_N"+to_string(i)); // Pulse amplitude (pedestal subtracted, ADC units)
		for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("energy_S"+to_string(i)); // Pulse amplitude (pedestal subtracted, ADC units)
	}
	
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("amp_N"+to_string(i)); // Pulse amplitude (pedestal subtracted, ADC units)
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("amp_S"+to_string(i)); // Pulse amplitude (pedestal subtracted, ADC units)
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("ped_N"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("ped_S"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("t_N"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("t_S"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("tcoarse_N"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("tcoarse_S"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("tfine_N"+to_string(i));// Number of pulses in event for this channel 
	for(int i=0; i<NChannels; ++i) locBranchRegister.Register_Single<Float_t>("tfine_S"+to_string(i));// Number of pulses in event for this channel 
	
	if(isRawMode) {
		locBranchRegister.Register_Single<Int_t>("nRawArrSize");
		for(int i =0; i < NChannels; ++i) locBranchRegister.Register_FundamentalArray<uint16_t>("raw_N" + to_string(i), "nRawArrSize", nRawArrSize);
		for(int i =0; i < NChannels; ++i) locBranchRegister.Register_FundamentalArray<uint16_t>("raw_S" + to_string(i), "nRawArrSize", nRawArrSize);
	}
	
	// Histograms to check matching between PS and beam photon tagger
	/// In terms of DeltaEnergy, DeltaTime
	/// Large peak at 0,0 indicates that beam photon matches between PS and tagger
	/// Early data quality checks showed that events with poorer matching are just as well resolved. No reason to add cut. (tagger is just inefficient, PS response same)
    // hPSTAGM_tdiffVsEdiff = new TH2F("PSTAGM_tdiffVsEdiff","PS pair - TAGM: PS-TAGM time difference vs. PS-TAGM energy difference;E(PS) - E(TAGM) [GeV];PSC/TAGM time difference [ns]",200,-1.0,1.0,1000,-200.,200.);
    // hPSTAGH_tdiffVsEdiff = new TH2F("PSTAGH_tdiffVsEdiff","PS pair - TAGH: PS-TAGH time difference vs. PS-TAGH energy difference;E(PS) - E(TAGH) [GeV];PSC/TAGH time difference [ns]",200,-1.0,1.0,1000,-200.,200.);
	
	dTreeInterface->Create_Branches(locBranchRegister);

  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventProcessor_BBCAL_tree::brun(JEventLoop *eventLoop, int32_t runnumber){

  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventProcessor_BBCAL_tree::evnt(JEventLoop *loop, uint64_t eventnumber)
{

	Bool_t read_out_event=false; // Only save event if at least one channel has amplitude > 10 above pedestal

	// Set up Pair Spectrometer vars (whether or not PS data read in)
	int tilel = -100;
	int tiler = -100;
	double tilelE = -100.;
	double tilerE = -100.;
	vector<const DPSPair*>   ps_pairs;
	vector<const DPSCPair*>   psc_pairs;
	// Tagger data
    vector<const DTAGHHit*> taghhits;
    vector<const DTAGMHit*> tagmhits;

	// Skip SYNC events, no useful data there
	if(eventnumber==0) return NOERROR; 
	

	if(!isCosmicSetup) {
		// TRIGGER: check that the correct trigger bit are found in event (from Sasha, logic not entirely clear to me)
		vector<const DL1Trigger *> l1trig;
		loop->Get(l1trig);
		unsigned int trig_bit[33];
		memset(trig_bit,0,sizeof(trig_bit));
		if( l1trig.size() > 0) {
			for(unsigned int bit = 0; bit < 32; bit++) trig_bit[bit+1] = (l1trig[0]->trig_mask & (1 << bit))    ? 1 : 0;       
		}
		// Removes most (or all(?)) with no PS pairs. Some passing cut still have no PS pair though.
		if(trig_bit[4] != 1) return NOERROR;
		
		// Retrieve PS data
		loop->Get(ps_pairs);
		loop->Get(psc_pairs);
		loop->Get(taghhits);
		loop->Get(tagmhits);
		
		// Looks like about about 18% here have 0 hits, 2% have 2+ hits.
		if(ps_pairs.size()!=1 || psc_pairs.size()!=1) {
			//if(print_verbose_out) cout << endl;
			return NOERROR;
		}
		tilel = ps_pairs[0]->ee.first->column; // left=north
		tiler = ps_pairs[0]->ee.second->column;// right=south
		tilelE = ps_pairs[0]->ee.first->E;
		tilerE = ps_pairs[0]->ee.second->E;
	}

	// Create arrays to store f250 readout data
	Float_t amp_N[NChannels], amp_S[NChannels]; 
	Float_t int_N[NChannels], int_S[NChannels]; 
	Float_t t_N[NChannels],   t_S[NChannels]; 
	uint32_t tcoarse_N[NChannels],   tcoarse_S[NChannels]; 
	uint32_t tfine_N[NChannels],   tfine_S[NChannels]; 
	Int_t nhit_N[NChannels],  nhit_S[NChannels]; 	
	Int_t ped_N[NChannels],   ped_S[NChannels]; 	
	memset(amp_N,0,sizeof(amp_N));   memset(amp_S,0,sizeof(amp_S)); 
	memset(int_N,0,sizeof(int_N));   memset(int_S,0,sizeof(int_S)); 
	memset(t_N,0,sizeof(t_N));       memset(t_S,0,sizeof(t_S)); 
	memset(tcoarse_N,0,sizeof(t_N)); memset(tcoarse_S,0,sizeof(t_S)); 
	memset(tfine_N,0,sizeof(t_N));   memset(tfine_S,0,sizeof(t_S)); 
	memset(nhit_N,0,sizeof(nhit_N)); memset(nhit_S,0,sizeof(nhit_S)); 
	memset(ped_N,0,sizeof(ped_N));   memset(ped_S,0,sizeof(ped_S)); 

	// Read in FADC data
	vector<const Df250PulseData*> pulse_data; 
	loop->Get(pulse_data);
	// verbose_counter++; // Only count up if we get past first two cuts

	for(unsigned int i = 0; i < pulse_data.size(); i++){

		int crate       =  pulse_data[i]->rocid;
		int slot        =  pulse_data[i]->slot;
		int channel     =  pulse_data[i]->channel;
		
		pair<int,int>  index_pair = getBBCALIndex(crate,slot,channel);
		int isSouth  = index_pair.first;
		int ch_index = index_pair.second;
		if(ch_index==-1) continue; // fADC250 pulse doesn't correspond to test stand

		// Calculate pedestal subtracted values
		uint32_t  coarse_time = pulse_data[i]->course_time; // Pretty sure "course" is a typo in the code repo...
		uint32_t  fine_time   = pulse_data[i]->fine_time;
		float pulse_time      = float(coarse_time)*4 + float(fine_time)*0.0625;
		int pulse_int         = pulse_data[i]->integral;
		int pulse_peak        = pulse_data[i]->pulse_peak;
		int pedestal          = pulse_data[i]->pedestal;
		float ped_sub_peak    = pulse_peak - (pedestal/4.);
		float ped_sub_int     = pulse_int  - (pedestal/4.)*(pulse_data[i]->nsamples_integral);
		
		if(ped_sub_peak<0.) ped_sub_peak=0;
		if(ped_sub_int<0.)  ped_sub_int=0;
		 
		// North channels: store fADC readouts
		if(isSouth==0) {
			nhit_N[ch_index]++;
			// Check for 2+ fADC250 pulses for this channel and event, skip all but first pulse
			// if(nhit_N[ch_index]!=1) {
				// Ignore second pulses
				// continue;
			// }
			// Fill output array with first fADC250 pulse encountered for this channel
			// else{
				if(ped_sub_peak>10.)  read_out_event=true;
				amp_N[ch_index]     = ped_sub_peak;
				int_N[ch_index]     = ped_sub_int;
				t_N[ch_index]       = pulse_time;
				tcoarse_N[ch_index] = coarse_time;
				tfine_N[ch_index]  = fine_time;
				ped_N[ch_index] = pedestal;
			// }
		}
		
		// South channels
		if(isSouth==1) {
			nhit_S[ch_index]++;
			// Check for 2+ fADC250 pulses for this channel and event, skip all but first pulse
			// if(nhit_S[ch_index]!=1) {
				// continue;
			// }
			// Fill output array with first fADC250 pulse encountered for this channel
			// else{
				if(ped_sub_peak>10.)  read_out_event=true;
				amp_S[ch_index] = ped_sub_peak;
				int_S[ch_index] = ped_sub_int;
				t_S[ch_index]       = pulse_time;
				tcoarse_S[ch_index] = coarse_time;
				tfine_S[ch_index]  = fine_time;
				ped_S[ch_index] = pedestal;
			// }
		}
	} // end of loop through Df250PulseData objects
	
	
	// Loop over raw waveforms, if data has them
	Int_t nhitRaw_N[NChannels], nhitRaw_S[NChannels];
	int16_t maxAmpRaw_N[NChannels], maxAmpRaw_S[NChannels];
	int16_t raw_N[NChannels][nRawArrSize], raw_S[NChannels][nRawArrSize];
	
	memset(maxAmpRaw_N,0,sizeof(maxAmpRaw_N));
	memset(maxAmpRaw_S,0,sizeof(maxAmpRaw_S));
	
	for (int ch = 0; ch < NChannels; ch++)
	{
		nhitRaw_N[ch]=0; nhitRaw_S[ch]=0;  
		for (int s = 0; s < nRawArrSize; s++){raw_N[ch][s] = 0; raw_S[ch][s] = 0;} 
	}
	
	
	
	
	if(isRawMode){
		vector<const Df250WindowRawData*> window_raw_data;
		loop->Get(window_raw_data);
				
		// cout << "number of raw mode objects this event: " << window_raw_data.size() << endl;
		for (unsigned int i = 0; i < window_raw_data.size(); i++)
		{
			int crate = window_raw_data[i]->rocid;
			int slot = window_raw_data[i]->slot;
			int channel = window_raw_data[i]->channel;
			
			
			// cout << "Rawmode crate: " << crate << " slot " << slot << " channel " << channel << endl;
			
			if(window_raw_data[i]->invalid_samples) continue;

			pair<int, int> index_pair = getBBCALIndex(crate,slot,channel);
			int isSouth  = index_pair.first;
			int ch_index = index_pair.second;
			if(ch_index==-1) continue;

			// cout << "CHANNEL BELONGS TO BABY BCAL! "<< endl;
			
			const vector<uint16_t> &samplesvector = window_raw_data[i]->samples;
			const int nsamples = samplesvector.size();
			
			// cout << "nsamples: " << samplesvector.size() << endl;
			
			if(isSouth) {
				nhitRaw_S[ch_index]++;
				auto max = *max_element(std::begin(samplesvector), std::end(samplesvector)); 
				if(maxAmpRaw_S[ch_index]<max) {
					maxAmpRaw_S[ch_index] = max;
					read_out_event=true;
					amp_S[ch_index]=max-ped_S[ch_index]/4.;
					for (int ii = 0; ii < nsamples; ii++) if(samplesvector[ii]>0) raw_S[ch_index][ii] = samplesvector[ii];
				}
			}
			
			if(!isSouth) {
				nhitRaw_N[ch_index]++;
				auto max = *max_element(std::begin(samplesvector), std::end(samplesvector)); 
				if(maxAmpRaw_N[ch_index]<max) {
					maxAmpRaw_N[ch_index] = max;
					read_out_event=true;
					amp_N[ch_index]=max-ped_N[ch_index]/4.;
					for (int ii = 0; ii < nsamples; ii++) if(samplesvector[ii]>0) raw_N[ch_index][ii] = samplesvector[ii];
				}
			}
			
		}
	}
		
	
	// Skip saving event if sum of amplitudes over channels is still essentially zero.
	if(!read_out_event) return NOERROR;
		
	japp->RootWriteLock();

	//FILL BRANCHES, THEN TTREE
	dTreeFillData.Fill_Single<int32_t>("run",loop->GetJEvent().GetRunNumber());
	dTreeFillData.Fill_Single<ULong64_t>("event",eventnumber);
	Float_t eshow_N=0.;    Float_t eshow_S=0.; 
	Float_t eshow_N_2S=0.; Float_t eshow_S_2S=0.; 
	for(int i=0; i<NChannels; ++i) {
		string i_str=to_string(i);
		// cout << "NEW AMP S IS: " << i << " " << amp_S[i] << endl;
		dTreeFillData.Fill_Single<Float_t>("amp_N"+i_str,     amp_N[i]);        dTreeFillData.Fill_Single<Float_t>("amp_S"+i_str,  amp_S[i]);
		// dTreeFillData.Fill_Single<Float_t>("int_N"+i_str,  int_N[i]);        dTreeFillData.Fill_Single<Float_t>("int_S"+i_str,  int_S[i]);
		dTreeFillData.Fill_Single<Float_t>("ped_N"+i_str,     ped_N[i]);        dTreeFillData.Fill_Single<Float_t>("ped_S"+i_str,  ped_S[i]);
		dTreeFillData.Fill_Single<Float_t>("t_N"+i_str,       t_N[i]);          dTreeFillData.Fill_Single<Float_t>("t_S"+i_str,    t_S[i]);
		dTreeFillData.Fill_Single<Float_t>("tcoarse_N"+i_str, tcoarse_N[i]);    dTreeFillData.Fill_Single<Float_t>("tcoarse_S"+i_str,    tcoarse_S[i]);
		dTreeFillData.Fill_Single<Float_t>("tfine_N"+i_str,   tfine_N[i]);      dTreeFillData.Fill_Single<Float_t>("tfine_S"+i_str,    tfine_S[i]);
		if(isRawMode) {
			dTreeFillData.Fill_Single<Int_t>("nRawArrSize", nRawArrSize);
			for (int ii = 0; ii < nRawArrSize; ii++)dTreeFillData.Fill_Array<uint16_t>("raw_N"+i_str, raw_N[i][ii], ii);
			for (int ii = 0; ii < nRawArrSize; ii++)dTreeFillData.Fill_Array<uint16_t>("raw_S"+i_str, raw_S[i][ii], ii);
		}
		
		if(!isRawMode) {
			eshow_N+=int_N[i]*C_N[i];
			eshow_S+=int_S[i]*C_S[i];
			if(int_N[i]>50. && int_S[i]>50.) {
				eshow_N_2S+=int_N[i]*C_N[i];
				eshow_S_2S+=int_S[i]*C_S[i];
			}
		}
		if(isRawMode) {
			eshow_N+=amp_N[i]*C_N[i];
			eshow_S+=amp_S[i]*C_S[i];
		}
		
	}


	if(myfile!=NULL) {
		for(int i=0; i<NChannels; ++i) {
			string i_str=to_string(i);
			if(!isRawMode) dTreeFillData.Fill_Single<Float_t>("energy_N"+i_str,  int_N[i]*C_N[i]); // This is actually energy/2, but when we sum over north+south it works out fine
			if(!isRawMode) dTreeFillData.Fill_Single<Float_t>("energy_S"+i_str,  int_S[i]*C_S[i]); // This is actually energy/2, but when we sum over north+south it works out fine
			if(isRawMode) dTreeFillData.Fill_Single<Float_t>("energy_N"+i_str,   amp_N[i]*C_N[i]*2.); // Different convention for cosmic data, true energy
			if(isRawMode) dTreeFillData.Fill_Single<Float_t>("energy_S"+i_str,   amp_S[i]*C_S[i]*2.); // Different convention for cosmic data, true energy
		}
		dTreeFillData.Fill_Single<Float_t>("energy_miniBCAL", eshow_N+eshow_S); // south Energy
		dTreeFillData.Fill_Single<Float_t>("energy_miniBCAL_2Sided", eshow_N_2S+eshow_S_2S); // south Energy
		dTreeFillData.Fill_Single<Float_t>("energy_miniBCAL_N", eshow_N); // south Energy
		dTreeFillData.Fill_Single<Float_t>("energy_miniBCAL_S", eshow_S); // south Energy
		dTreeFillData.Fill_Single<Float_t>("energy_miniBCAL_N_2Sided", eshow_N_2S); // south Energy
		dTreeFillData.Fill_Single<Float_t>("energy_miniBCAL_S_2Sided", eshow_S_2S); // south Energy
	}
	
	if(!isCosmicSetup) {
		dTreeFillData.Fill_Single<Int_t>("PS_NumPairs",ps_pairs.size());
		dTreeFillData.Fill_Single<Int_t>("PSC_NumPairs",psc_pairs.size());
		dTreeFillData.Fill_Single<Int_t>("PS_column_left",tilel); // left=north
		dTreeFillData.Fill_Single<Int_t>("PS_column_right",tiler);// right=south
		dTreeFillData.Fill_Single<Double_t>("PS_column_right_E", tilerE); // south Energy
		dTreeFillData.Fill_Single<Double_t>("PS_column_left_E",  tilelE); // south Energy
		
		for(int i=0; i<NChannels; ++i) {
			string i_str=to_string(i);
			dTreeFillData.Fill_Single<Float_t>("int_N"+i_str, int_N[i]); 
			dTreeFillData.Fill_Single<Float_t>("int_S"+i_str, int_S[i]); 
		}
	}
		
/* 		// TAGH/TAGM for energy+time matching
		// Jon was concerned about potential junk hits in PS, so looked into requiring that different subsystem show matching energy sum
		// In practice, this cut wasn't needed, didn't change PS quality at all (checks for good PS hits already sufficient).
		Double_t  max_E_diff = 0.2;
		Int_t NTAGH_PASS_ECUT=0;
		Int_t NTAGM_PASS_ECUT=0;
		const DPSCHit* clhit = psc_pairs[0]->ee.first; // left hit in coarse PS
		// Loop over TAGM hits
		for (unsigned int i=0; i < tagmhits.size(); i++) {
			const DTAGMHit* tag = tagmhits[i];
			if (!tag->has_TDC||!tag->has_fADC) continue;
			if (tag->row!=0) continue; // Copypasta from PSPair_online plugin... Don't actually know what this does
			Double_t Ediff = tilelE+tilerE-tag->E;
			Double_t tdiff = clhit->t-tag->t;
			if(fabs(Ediff)<max_E_diff) {
				dTreeFillData.Fill_Array<Double_t>("Ediff_tagm", Ediff, NTAGM_PASS_ECUT);
				dTreeFillData.Fill_Array<Double_t>("tdiff_tagm", tdiff, NTAGM_PASS_ECUT);
				NTAGM_PASS_ECUT++;
			}
			hPSTAGM_tdiffVsEdiff->Fill(Ediff,tdiff);
		}
		// dTreeFillData.Fill_Single<Int_t>("NTAGM_PASSING",  NTAGM_PASS_ECUT);
 		
		
		// Loop over TAGH hits
		for (unsigned int i=0; i < taghhits.size(); i++) {
			const DTAGHHit* tag = taghhits[i];
			if (!tag->has_TDC||!tag->has_fADC) continue;
			Double_t Ediff = tilelE+tilerE-tag->E;
			Double_t tdiff = clhit->t-tag->t;
			if(fabs(Ediff)<max_E_diff) {
				dTreeFillData.Fill_Array<Double_t>("Ediff_tagh", Ediff, NTAGH_PASS_ECUT);
				dTreeFillData.Fill_Array<Double_t>("tdiff_tagh", tdiff, NTAGH_PASS_ECUT);
				NTAGH_PASS_ECUT++;
			}
			hPSTAGH_tdiffVsEdiff->Fill(Ediff,tdiff);
		}
		dTreeFillData.Fill_Single<Int_t>("NTAGH_PASSING",  NTAGH_PASS_ECUT);
	}
*/
		
	dTreeInterface->Fill(dTreeFillData);
	japp->RootUnLock();


	return NOERROR;
	
}

//------------------
// erun
//------------------
jerror_t DEventProcessor_BBCAL_tree::erun(void)
{
	
	// Any final calculations on histograms (like dividing them)
	// should be done here. This may get called more than once.
	return NOERROR;
}

// void DEventProcessor_BBCAL_tree::ParseGainFile(ifstream input) {
// void DEventProcessor_BBCAL_tree::ParseGainFile(const char* fname) {
void DEventProcessor_BBCAL_tree::ParseGainFile(string fname) {
	
 	// Copypasta to make a 2D vector of lines/rows (assuming space separator)
	typedef vector<vector<string> > DRows;
	vector<vector<string> > parsed;
	ifstream input(fname);
	char const row_delim = '\n';
	char const field_delim = ' ';
	for (string row; getline(input, row, row_delim); ) {
	  parsed.push_back(DRows::value_type());
	  istringstream ss(row);
	  for (string field; getline(ss, field, field_delim); ) {
		parsed.back().push_back(field.c_str());
	  }
	}
	input.close();
	
	// Loop over N+S channel
	for(size_t i=0; i<16; i++) {
		string northname="N"+to_string(i)+":";
		string southname="S"+to_string(i)+":";
		for(size_t ii=0; ii<parsed.size(); ii++) {
			if(parsed[ii].size()!=0) {
				string first_in_line=parsed[ii][0];
				if(first_in_line==northname) C_N[i]=atof(parsed[ii][1].c_str());
				if(first_in_line==southname) C_S[i]=atof(parsed[ii][1].c_str());
			}
		}
	}
	
return;
}



//------------------
// fini
//------------------
jerror_t DEventProcessor_BBCAL_tree::fini(void)
{
	if(dTreeInterface != NULL) delete dTreeInterface; // SAVES FILES, TREES during deletion step
	
	return NOERROR;
}

