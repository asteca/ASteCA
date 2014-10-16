# isochrones folder

At least one folder with Girardi isochrone files should be located here.

**ASteCA** supports any set of [CMD Girardi et al. theoretical isochrones](http://stev.oapd.inaf.it/cgi-bin/cmd).

The convention for naming the folders is to use the identificators defined in the CMD site for each photometric system:

* `2mass_spitzer`:  2MASS + Spitzer (IRAC+MIPS)
* `2mass_spitzer_wise`:  2MASS + Spitzer (IRAC+MIPS) + WISE
* `2mass`:  2MASS JHK_s
* `ogle_2mass_spitzer`:  OGLE + 2MASS + Spitzer (IRAC+MIPS)
* `2mass_spitzer_wise_washington_ddo51`: 2MASS+Spitzer+WISE+Washington+DDO51
* `ubvrijhk`: UBVRIJHK (cf. Maiz-Apellaniz 2006 + Bessell 1990)
* `bessell`: UBVRIJHKLMN (cf. Bessell 1990 + Bessell & Brett 1988)
* `akari`: AKARI
* `batc`: BATC
* `megacam_wircam`: CFHT Megacam + Wircam (all ABmags)
* `wircam`: CFHT Wircam
* `megacam`: CFHT/Megacam u*g'r'i'z'
* `ciber`: CIBER
* `dcmc`: DCMC
* `decam`: DECAM (ABmags)
* `decam_vista`: DECAM ugrizY (ABmags) + VISTA ZYJHK_s (Vegamags)
* `denis`: DENIS
* `dmc14`: DMC 14 filters
* `dmc15`: DMC 15 filters
* `eis`: ESO/EIS (WFI UBVRIZ + SOFI JHK)
* `wfi`: ESO/WFI
* `wfi_sofi`: ESO/WFI+SOFI
* `wfi2`: ESO/WFI2
* `galex`: GALEX FUV+NUV (Vegamag) + Johnson's UBV
* `galex_sloan`: GALEX FUV+NUV + SDSS ugriz (all ABmags) 
* `UVbright`: HST+GALEX+Swift/UVOT UV filters
* `acs_hrc`: HST/ACS HRC
* `acs_wfc`: HST/ACS WFC
* `nicmosab`: HST/NICMOS AB
* `nicmosst`: HST/NICMOS ST
* `nicmosvega`: HST/NICMOS vega
* `stis`: HST/STIS imaging mode
* `wfc3ir`: HST/WFC3 IR channel (final throughputs)
* `wfc3uvis1`: HST/WFC3 UVIS channel, chip 1 (final throughputs)
* `wfc3uvis2`: HST/WFC3 UVIS channel, chip 2 (final throughputs)
* `wfc3_verywide`: HST/WFC3 long-pass and extremely wide filters (UVIS1, final throughputs)
* `wfc3_medium`: HST/WFC3 medium filters (UVIS1+IR, final throughputs)
* `wfc3_wide`: HST/WFC3 wide filters (UVIS1+IR, final throughputs)
* `wfpc2`: HST/WFPC2 (Vegamag, cf. Holtzman et al. 1995)
* `int_wfc`: INT/WFC (Vegamag)
* `iphas`: IPHAS
* `kepler`: Kepler + SDSS griz + DDO51 (in ABmags)
* `kepler_2mass`: Kepler + SDSS griz + DDO51 (in ABmags) + 2MASS (~Vegamag)
* `lbt_lbc`: LBT/LBC (Vegamag)
* `noao_ctio_mosaic2`: NOAO/CTIO/MOSAIC2 (Vegamag)
* `ogle`: OGLE-II
* `panstarrs1`: Pan-STARRS1
* `sloan`: SDSS ugriz
* `sloan_2mass`: SDSS ugriz + 2MASS JHK_s
* `sloan_ukidss`: SDSS ugriz + UKIDSS ZYJHK
* `swift_uvot`: SWIFT/UVOT UVW2, UVM2, UVW1,u (Vegamag) 
* `spitzer`: Spitzer IRAC+MIPS
* `stroemgren`: Stroemgren-Crawford
* `suprimecam`: Subaru/Suprime-Cam (ABmags)
* `TESS_2mass_kepler`: TESS + 2MASS (Vegamags) + Kepler + SDSS griz + DDO51 (in ABmags)
* `tycho2`: Tycho V_TB_T
* `ukidss`: UKIDSS ZYJHK (Vegamag)
* `visir`: VISIR
* `vista`: VISTA ZYJHK_s (Vegamag)
* `vst_omegacam`: VST/OMEGACAM (Vegamag)
* `vilnius`: Vilnius
* `washington`: Washington CMT_1T_2
* `washington_ddo51`: Washington CMT_1T_2 + DDO51

This file can be safely deleted.