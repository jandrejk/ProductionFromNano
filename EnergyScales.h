const struct ES{
	struct Electron{
		float oneProng0p0 = 0.0;
		float oneProng1p0 = 0.0;
		float threeProng0p0 = 0.0;
		float uncertaintyShift1p   = 0.02;
		float uncertaintyShift1p1p = 0.02;
		float uncertaintyShift3p   = 0.02;
		float uncertaintyBarrel    = 0.02;
		float uncertaintyEndcap    = 0.02;
	} Electron;

	struct Muon{
		float oneProng0p0 = 0.0;
		float oneProng1p0 = 0.0;
		float threeProng0p0 = 0.0;
		float uncertaintyShift1p   = 0.02;
		float uncertaintyShift1p1p = 0.02;
		float uncertaintyShift3p   = 0.02;
	} Muon;
	// taken from here https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
	struct Tau{
		float oneProng0p0 =   -0.013;
		float oneProng1p0 =   -0.005;
		float threeProng0p0 = -0.012;
		float uncertaintyShift1p   = 0.011;
		float uncertaintyShift1p1p = 0.009;
		float uncertaintyShift3p   = 0.008;
	} Tau;
} ES;