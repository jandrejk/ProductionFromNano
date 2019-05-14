namespace LeptonCuts
{

	// taken from https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2016#Baseline_mu_tau_h
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	const struct Baseline
	{
		struct Electron
		{
			float pt  = 25.0;
			float eta = 2.1;
		}Electron;

		struct Muon
		{
			float pt  = 20.0;
			float eta = 2.1;
		}Muon;

		struct Tau
		{
			int Additional = 8;
			struct SemiLep
			{
				float pt  = 20.0;
				float eta = 2.3;
				int bitmask = 16;
			}SemiLep;

			struct FullHad
			{
				float lead_pt  = 40.0;
				float sublead_pt  = 40.0;
				float eta = 2.1;
				int lead_bitmask = 32;
				int sub_bitmask  = 64;
				
			}FullHad;
		}Tau;

		int bitmask = 1;

	}Baseline;
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	const struct Di
	{
		struct Electron
		{
			float pt  = 15.0;
			float eta = 2.4;
		}Electron;

		struct Muon
		{
			float pt  = 15.0;
			float eta = 2.4;
		}Muon;

		int bitmask = 2;

	}Di;
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	const struct Extra //Third lepton vetoes
	{
		struct Electron
		{
			float pt  = 10.0;
			float eta = 2.5;
		}Electron;

		struct Muon
		{
			float pt  = 10.0;
			float eta = 2.4;
		}Muon;

		int bitmask = 4;

	}Extra;
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	const struct Additional
	{
		struct Electron
		{
			float pt  = 13.0;
			float eta = 2.1;
		}Electron;

		struct Muon
		{
			float pt  = 5.0;
			float eta = 2.1;
		}Muon;

		int bitmask = 8;

	}Additional;
	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
}
