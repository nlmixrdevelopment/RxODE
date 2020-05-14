//#undef NDEBUG
#include <R.h>
#include <stan/math.hpp>
#include "../inst/include/RxODE.h"
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("RxODE", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif
extern "C" double linCmtB(rx_solve *rx, unsigned int id, double t, int linCmt,
			  int i_cmt, int trans, int val,
			  double dd_p1, double dd_v1,
			  double dd_p2, double dd_p3,
			  double dd_p4, double dd_p5,
			  double dd_ka,
			  double dd_tlag, double dd_tlag2,
			  double dd_F, double dd_F2,
			  double dd_rate, double dd_dur,
			  double dd_rate2, double dd_dur2){
  return 0.0;
}
