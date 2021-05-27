#ifndef __cbindThetaOmega_H__
#define __cbindThetaOmega_H__
#if defined(__cplusplus)
List cbindThetaOmega(RObject inputParametersRO, List& individualParameters);
extern "C" {
#endif

SEXP _rxCbindStudyIndividual(SEXP inputParameters, SEXP individualParameters);
#if defined(__cplusplus)
}
#endif
#endif
