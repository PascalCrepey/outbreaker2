// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_are_possible_ancestors
std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i);
RcppExport SEXP outbreaker2_cpp_are_possible_ancestors(SEXP t_infSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type t_inf(t_infSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_are_possible_ancestors(t_inf, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_sample1
size_t cpp_sample1(std::vector<int> x);
RcppExport SEXP outbreaker2_cpp_sample1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sample1(x));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pick_possible_ancestor
size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i);
RcppExport SEXP outbreaker2_cpp_pick_possible_ancestor(SEXP t_infSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type t_inf(t_infSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pick_possible_ancestor(t_inf, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_find_descendents
Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i);
RcppExport SEXP outbreaker2_cpp_find_descendents(SEXP alphaSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_find_descendents(alpha, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_find_local_cases
Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i);
RcppExport SEXP outbreaker2_cpp_find_local_cases(SEXP alphaSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_find_local_cases(alpha, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_swap_cases
Rcpp::List cpp_swap_cases(Rcpp::List param, size_t i);
RcppExport SEXP outbreaker2_cpp_swap_cases(SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_swap_cases(param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_genetic
double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i);
RcppExport SEXP outbreaker2_cpp_ll_genetic(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_genetic(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_timing_infections
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i);
RcppExport SEXP outbreaker2_cpp_ll_timing_infections(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_timing_infections(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_timing_sampling
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i);
RcppExport SEXP outbreaker2_cpp_ll_timing_sampling(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_timing_sampling(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_reporting
double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i);
RcppExport SEXP outbreaker2_cpp_ll_reporting(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_reporting(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_timing
double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i);
RcppExport SEXP outbreaker2_cpp_ll_timing(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_timing(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_all
double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i);
RcppExport SEXP outbreaker2_cpp_ll_all(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_all(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_mu
Rcpp::List cpp_move_mu(Rcpp::List data, Rcpp::List param, Rcpp::List config);
RcppExport SEXP outbreaker2_cpp_move_mu(SEXP dataSEXP, SEXP paramSEXP, SEXP configSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_mu(data, param, config));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_t_inf
Rcpp::List cpp_move_t_inf(Rcpp::List data, Rcpp::List param);
RcppExport SEXP outbreaker2_cpp_move_t_inf(SEXP dataSEXP, SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_t_inf(data, param));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_alpha
Rcpp::List cpp_move_alpha(Rcpp::List data, Rcpp::List param);
RcppExport SEXP outbreaker2_cpp_move_alpha(SEXP dataSEXP, SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_alpha(data, param));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_swap_cases
Rcpp::List cpp_move_swap_cases(Rcpp::List data, Rcpp::List param);
RcppExport SEXP outbreaker2_cpp_move_swap_cases(SEXP dataSEXP, SEXP paramSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_swap_cases(data, param));
    return rcpp_result_gen;
END_RCPP
}
