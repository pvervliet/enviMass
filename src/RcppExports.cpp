// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getEIC_new
SEXP getEIC_new(SEXP mz, SEXP RT, SEXP intens, SEXP orderedint, SEXP orderedret, SEXP dmzdens, SEXP ppm2, SEXP drtdens, SEXP merged2);
RcppExport SEXP _enviMass_getEIC_new(SEXP mzSEXP, SEXP RTSEXP, SEXP intensSEXP, SEXP orderedintSEXP, SEXP orderedretSEXP, SEXP dmzdensSEXP, SEXP ppm2SEXP, SEXP drtdensSEXP, SEXP merged2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mz(mzSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intens(intensSEXP);
    Rcpp::traits::input_parameter< SEXP >::type orderedint(orderedintSEXP);
    Rcpp::traits::input_parameter< SEXP >::type orderedret(orderedretSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dmzdens(dmzdensSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ppm2(ppm2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type drtdens(drtdensSEXP);
    Rcpp::traits::input_parameter< SEXP >::type merged2(merged2SEXP);
    rcpp_result_gen = Rcpp::wrap(getEIC_new(mz, RT, intens, orderedint, orderedret, dmzdens, ppm2, drtdens, merged2));
    return rcpp_result_gen;
END_RCPP
}
// series_relat
SEXP series_relat(SEXP homol_peaks_relat, SEXP range_mz, SEXP range_RT);
RcppExport SEXP _enviMass_series_relat(SEXP homol_peaks_relatSEXP, SEXP range_mzSEXP, SEXP range_RTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type homol_peaks_relat(homol_peaks_relatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type range_mz(range_mzSEXP);
    Rcpp::traits::input_parameter< SEXP >::type range_RT(range_RTSEXP);
    rcpp_result_gen = Rcpp::wrap(series_relat(homol_peaks_relat, range_mz, range_RT));
    return rcpp_result_gen;
END_RCPP
}
// moving_count
SEXP moving_count(SEXP homol_peaks_relat, SEXP deldel);
RcppExport SEXP _enviMass_moving_count(SEXP homol_peaks_relatSEXP, SEXP deldelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type homol_peaks_relat(homol_peaks_relatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type deldel(deldelSEXP);
    rcpp_result_gen = Rcpp::wrap(moving_count(homol_peaks_relat, deldel));
    return rcpp_result_gen;
END_RCPP
}
// compare
SEXP compare(SEXP a, SEXP b, SEXP results);
RcppExport SEXP _enviMass_compare(SEXP aSEXP, SEXP bSEXP, SEXP resultsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type b(bSEXP);
    Rcpp::traits::input_parameter< SEXP >::type results(resultsSEXP);
    rcpp_result_gen = Rcpp::wrap(compare(a, b, results));
    return rcpp_result_gen;
END_RCPP
}
// correct_intens
SEXP correct_intens(SEXP corfac, SEXP sampleID, SEXP intens, SEXP sampleID_peak);
RcppExport SEXP _enviMass_correct_intens(SEXP corfacSEXP, SEXP sampleIDSEXP, SEXP intensSEXP, SEXP sampleID_peakSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type corfac(corfacSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sampleID(sampleIDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intens(intensSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sampleID_peak(sampleID_peakSEXP);
    rcpp_result_gen = Rcpp::wrap(correct_intens(corfac, sampleID, intens, sampleID_peak));
    return rcpp_result_gen;
END_RCPP
}
// metagroup
SEXP metagroup(SEXP proffrom, /* must be sorted */                             SEXP profto);
RcppExport SEXP _enviMass_metagroup(SEXP proffromSEXP, SEXP proftoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type proffrom(proffromSEXP);
    Rcpp::traits::input_parameter< /* must be sorted */                             SEXP >::type profto(proftoSEXP);
    rcpp_result_gen = Rcpp::wrap(metagroup(proffrom, profto));
    return rcpp_result_gen;
END_RCPP
}
// profpeakprof
SEXP profpeakprof(SEXP ProPeak_pro, SEXP ProPeak_peak, SEXP Peak_peak1, SEXP Peak_peak2, SEXP Peak_score, SEXP PeakPro);
RcppExport SEXP _enviMass_profpeakprof(SEXP ProPeak_proSEXP, SEXP ProPeak_peakSEXP, SEXP Peak_peak1SEXP, SEXP Peak_peak2SEXP, SEXP Peak_scoreSEXP, SEXP PeakProSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type ProPeak_pro(ProPeak_proSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ProPeak_peak(ProPeak_peakSEXP);
    Rcpp::traits::input_parameter< SEXP >::type Peak_peak1(Peak_peak1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Peak_peak2(Peak_peak2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type Peak_score(Peak_scoreSEXP);
    Rcpp::traits::input_parameter< SEXP >::type PeakPro(PeakProSEXP);
    rcpp_result_gen = Rcpp::wrap(profpeakprof(ProPeak_pro, ProPeak_peak, Peak_peak1, Peak_peak2, Peak_score, PeakPro));
    return rcpp_result_gen;
END_RCPP
}
// mergeProfiles
SEXP mergeProfiles(SEXP mz_lower, SEXP mz_upper, SEXP RT_lower, SEXP RT_upper, SEXP intens, SEXP sam, SEXP orderedint, SEXP orderedsam, SEXP supress);
RcppExport SEXP _enviMass_mergeProfiles(SEXP mz_lowerSEXP, SEXP mz_upperSEXP, SEXP RT_lowerSEXP, SEXP RT_upperSEXP, SEXP intensSEXP, SEXP samSEXP, SEXP orderedintSEXP, SEXP orderedsamSEXP, SEXP supressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mz_lower(mz_lowerSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mz_upper(mz_upperSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RT_lower(RT_lowerSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RT_upper(RT_upperSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intens(intensSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sam(samSEXP);
    Rcpp::traits::input_parameter< SEXP >::type orderedint(orderedintSEXP);
    Rcpp::traits::input_parameter< SEXP >::type orderedsam(orderedsamSEXP);
    Rcpp::traits::input_parameter< SEXP >::type supress(supressSEXP);
    rcpp_result_gen = Rcpp::wrap(mergeProfiles(mz_lower, mz_upper, RT_lower, RT_upper, intens, sam, orderedint, orderedsam, supress));
    return rcpp_result_gen;
END_RCPP
}
// neighbour
SEXP neighbour(SEXP mz, /* must be sorted */                        SEXP rt, SEXP sample, SEXP maxsample, SEXP ppm, SEXP dmz, SEXP drt);
RcppExport SEXP _enviMass_neighbour(SEXP mzSEXP, SEXP rtSEXP, SEXP sampleSEXP, SEXP maxsampleSEXP, SEXP ppmSEXP, SEXP dmzSEXP, SEXP drtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mz(mzSEXP);
    Rcpp::traits::input_parameter< /* must be sorted */                        SEXP >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type maxsample(maxsampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dmz(dmzSEXP);
    Rcpp::traits::input_parameter< SEXP >::type drt(drtSEXP);
    rcpp_result_gen = Rcpp::wrap(neighbour(mz, rt, sample, maxsample, ppm, dmz, drt));
    return rcpp_result_gen;
END_RCPP
}
// agglom
SEXP agglom(SEXP mz, /* must be sorted */                        SEXP rt, SEXP sample, SEXP ppm, SEXP dmz, SEXP drt);
RcppExport SEXP _enviMass_agglom(SEXP mzSEXP, SEXP rtSEXP, SEXP sampleSEXP, SEXP ppmSEXP, SEXP dmzSEXP, SEXP drtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mz(mzSEXP);
    Rcpp::traits::input_parameter< /* must be sorted */                        SEXP >::type rt(rtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sample(sampleSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dmz(dmzSEXP);
    Rcpp::traits::input_parameter< SEXP >::type drt(drtSEXP);
    rcpp_result_gen = Rcpp::wrap(agglom(mz, rt, sample, ppm, dmz, drt));
    return rcpp_result_gen;
END_RCPP
}
// indexed
SEXP indexed(SEXP index, /* must be sorted */                    SEXP maxindex, SEXP many);
RcppExport SEXP _enviMass_indexed(SEXP indexSEXP, SEXP maxindexSEXP, SEXP manySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type index(indexSEXP);
    Rcpp::traits::input_parameter< /* must be sorted */                    SEXP >::type maxindex(maxindexSEXP);
    Rcpp::traits::input_parameter< SEXP >::type many(manySEXP);
    rcpp_result_gen = Rcpp::wrap(indexed(index, maxindex, many));
    return rcpp_result_gen;
END_RCPP
}
// fill_timeset
SEXP fill_timeset(SEXP timeset, SEXP sampleID, SEXP intensity, SEXP lengtimeset);
RcppExport SEXP _enviMass_fill_timeset(SEXP timesetSEXP, SEXP sampleIDSEXP, SEXP intensitySEXP, SEXP lengtimesetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type timeset(timesetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sampleID(sampleIDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< SEXP >::type lengtimeset(lengtimesetSEXP);
    rcpp_result_gen = Rcpp::wrap(fill_timeset(timeset, sampleID, intensity, lengtimeset));
    return rcpp_result_gen;
END_RCPP
}
// meandel
SEXP meandel(SEXP timeset, SEXP subit, SEXP subrat, SEXP numtime, SEXP getwhat, SEXP lags, SEXP threshold, SEXP notrend);
RcppExport SEXP _enviMass_meandel(SEXP timesetSEXP, SEXP subitSEXP, SEXP subratSEXP, SEXP numtimeSEXP, SEXP getwhatSEXP, SEXP lagsSEXP, SEXP thresholdSEXP, SEXP notrendSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type timeset(timesetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subit(subitSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subrat(subratSEXP);
    Rcpp::traits::input_parameter< SEXP >::type numtime(numtimeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type getwhat(getwhatSEXP);
    Rcpp::traits::input_parameter< SEXP >::type lags(lagsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< SEXP >::type notrend(notrendSEXP);
    rcpp_result_gen = Rcpp::wrap(meandel(timeset, subit, subrat, numtime, getwhat, lags, threshold, notrend));
    return rcpp_result_gen;
END_RCPP
}
// intdiff
SEXP intdiff(SEXP timeset, SEXP subit, SEXP subrat, SEXP numtime, SEXP getwhat);
RcppExport SEXP _enviMass_intdiff(SEXP timesetSEXP, SEXP subitSEXP, SEXP subratSEXP, SEXP numtimeSEXP, SEXP getwhatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type timeset(timesetSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subit(subitSEXP);
    Rcpp::traits::input_parameter< SEXP >::type subrat(subratSEXP);
    Rcpp::traits::input_parameter< SEXP >::type numtime(numtimeSEXP);
    Rcpp::traits::input_parameter< SEXP >::type getwhat(getwhatSEXP);
    rcpp_result_gen = Rcpp::wrap(intdiff(timeset, subit, subrat, numtime, getwhat));
    return rcpp_result_gen;
END_RCPP
}
// plot_prof
SEXP plot_prof(SEXP RTlim_low, SEXP RTlim_up, SEXP mzlim_low, SEXP mzlim_up, SEXP mz, SEXP RT, SEXP intensity, SEXP sampleID, SEXP color1, SEXP color2, SEXP whatcolor);
RcppExport SEXP _enviMass_plot_prof(SEXP RTlim_lowSEXP, SEXP RTlim_upSEXP, SEXP mzlim_lowSEXP, SEXP mzlim_upSEXP, SEXP mzSEXP, SEXP RTSEXP, SEXP intensitySEXP, SEXP sampleIDSEXP, SEXP color1SEXP, SEXP color2SEXP, SEXP whatcolorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type RTlim_low(RTlim_lowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RTlim_up(RTlim_upSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mzlim_low(mzlim_lowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mzlim_up(mzlim_upSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mz(mzSEXP);
    Rcpp::traits::input_parameter< SEXP >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< SEXP >::type sampleID(sampleIDSEXP);
    Rcpp::traits::input_parameter< SEXP >::type color1(color1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type color2(color2SEXP);
    Rcpp::traits::input_parameter< SEXP >::type whatcolor(whatcolorSEXP);
    rcpp_result_gen = Rcpp::wrap(plot_prof(RTlim_low, RTlim_up, mzlim_low, mzlim_up, mz, RT, intensity, sampleID, color1, color2, whatcolor));
    return rcpp_result_gen;
END_RCPP
}
// binRT_prof
SEXP binRT_prof(SEXP RT, SEXP intensity, SEXP binRT, SEXP colorit, SEXP what);
RcppExport SEXP _enviMass_binRT_prof(SEXP RTSEXP, SEXP intensitySEXP, SEXP binRTSEXP, SEXP coloritSEXP, SEXP whatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type RT(RTSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< SEXP >::type binRT(binRTSEXP);
    Rcpp::traits::input_parameter< SEXP >::type colorit(coloritSEXP);
    Rcpp::traits::input_parameter< SEXP >::type what(whatSEXP);
    rcpp_result_gen = Rcpp::wrap(binRT_prof(RT, intensity, binRT, colorit, what));
    return rcpp_result_gen;
END_RCPP
}
// binmz_prof
SEXP binmz_prof(SEXP mz, SEXP intensity, SEXP binmzs, SEXP colorit);
RcppExport SEXP _enviMass_binmz_prof(SEXP mzSEXP, SEXP intensitySEXP, SEXP binmzsSEXP, SEXP coloritSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mz(mzSEXP);
    Rcpp::traits::input_parameter< SEXP >::type intensity(intensitySEXP);
    Rcpp::traits::input_parameter< SEXP >::type binmzs(binmzsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type colorit(coloritSEXP);
    rcpp_result_gen = Rcpp::wrap(binmz_prof(mz, intensity, binmzs, colorit));
    return rcpp_result_gen;
END_RCPP
}
// extractProfiles
IntegerVector extractProfiles(NumericMatrix peaks, IntegerVector in_order, double dmass, bool ppm, double dret);
RcppExport SEXP _enviMass_extractProfiles(SEXP peaksSEXP, SEXP in_orderSEXP, SEXP dmassSEXP, SEXP ppmSEXP, SEXP dretSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type peaks(peaksSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type in_order(in_orderSEXP);
    Rcpp::traits::input_parameter< double >::type dmass(dmassSEXP);
    Rcpp::traits::input_parameter< bool >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< double >::type dret(dretSEXP);
    rcpp_result_gen = Rcpp::wrap(extractProfiles(peaks, in_order, dmass, ppm, dret));
    return rcpp_result_gen;
END_RCPP
}
// extractProfiles_replicates
IntegerVector extractProfiles_replicates(NumericMatrix peaks, IntegerVector in_order, double dmass, bool ppm, double dret, IntegerVector pregroup);
RcppExport SEXP _enviMass_extractProfiles_replicates(SEXP peaksSEXP, SEXP in_orderSEXP, SEXP dmassSEXP, SEXP ppmSEXP, SEXP dretSEXP, SEXP pregroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type peaks(peaksSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type in_order(in_orderSEXP);
    Rcpp::traits::input_parameter< double >::type dmass(dmassSEXP);
    Rcpp::traits::input_parameter< bool >::type ppm(ppmSEXP);
    Rcpp::traits::input_parameter< double >::type dret(dretSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type pregroup(pregroupSEXP);
    rcpp_result_gen = Rcpp::wrap(extractProfiles_replicates(peaks, in_order, dmass, ppm, dret, pregroup));
    return rcpp_result_gen;
END_RCPP
}
// while_checked
List while_checked(List check_nodes, NumericMatrix pattern_compound, NumericMatrix peaks, double RT_tol_inside, double int_tol);
RcppExport SEXP _enviMass_while_checked(SEXP check_nodesSEXP, SEXP pattern_compoundSEXP, SEXP peaksSEXP, SEXP RT_tol_insideSEXP, SEXP int_tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type check_nodes(check_nodesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pattern_compound(pattern_compoundSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type peaks(peaksSEXP);
    Rcpp::traits::input_parameter< double >::type RT_tol_inside(RT_tol_insideSEXP);
    Rcpp::traits::input_parameter< double >::type int_tol(int_tolSEXP);
    rcpp_result_gen = Rcpp::wrap(while_checked(check_nodes, pattern_compound, peaks, RT_tol_inside, int_tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_enviMass_getEIC_new", (DL_FUNC) &_enviMass_getEIC_new, 9},
    {"_enviMass_series_relat", (DL_FUNC) &_enviMass_series_relat, 3},
    {"_enviMass_moving_count", (DL_FUNC) &_enviMass_moving_count, 2},
    {"_enviMass_compare", (DL_FUNC) &_enviMass_compare, 3},
    {"_enviMass_correct_intens", (DL_FUNC) &_enviMass_correct_intens, 4},
    {"_enviMass_metagroup", (DL_FUNC) &_enviMass_metagroup, 2},
    {"_enviMass_profpeakprof", (DL_FUNC) &_enviMass_profpeakprof, 6},
    {"_enviMass_mergeProfiles", (DL_FUNC) &_enviMass_mergeProfiles, 9},
    {"_enviMass_neighbour", (DL_FUNC) &_enviMass_neighbour, 7},
    {"_enviMass_agglom", (DL_FUNC) &_enviMass_agglom, 6},
    {"_enviMass_indexed", (DL_FUNC) &_enviMass_indexed, 3},
    {"_enviMass_fill_timeset", (DL_FUNC) &_enviMass_fill_timeset, 4},
    {"_enviMass_meandel", (DL_FUNC) &_enviMass_meandel, 8},
    {"_enviMass_intdiff", (DL_FUNC) &_enviMass_intdiff, 5},
    {"_enviMass_plot_prof", (DL_FUNC) &_enviMass_plot_prof, 11},
    {"_enviMass_binRT_prof", (DL_FUNC) &_enviMass_binRT_prof, 5},
    {"_enviMass_binmz_prof", (DL_FUNC) &_enviMass_binmz_prof, 4},
    {"_enviMass_extractProfiles", (DL_FUNC) &_enviMass_extractProfiles, 5},
    {"_enviMass_extractProfiles_replicates", (DL_FUNC) &_enviMass_extractProfiles_replicates, 6},
    {"_enviMass_while_checked", (DL_FUNC) &_enviMass_while_checked, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_enviMass(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
