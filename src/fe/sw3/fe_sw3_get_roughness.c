/* Given a depth and string id, returns a drag coefficient */
/* Link to the friction library */

#include "global_header.h"
#include "friction_ext.h"

double fe_sw3_get_roughness(
  SMODEL *mod,
  double depth,         /* local depth */
  int id           /* string id */
)
{
  double temp = 0.;
  double zero_drag = 0.0;
  double c_eda = 0.;

  if (mod->str_values[id].fterms.eqrheight > SMALL) {
    c_eda = fr_equiv_depavg_shstr_coef(depth, mod->str_values[id].fterms.eqrheight);
  }
  else {
    c_eda = 0.;
  }

  if (mod->str_values[id].fterms.sav_flag == YES) {
    temp = fr_sav_drag_coef(depth, mod->str_values[id].fterms.eqrheight, mod->str_values[id].fterms.hghtstem);
    return temp * c_eda;
  }
  else if (mod->str_values[id].fterms.urv_flag == YES) {
    temp = fr_urv_drag_coef(depth, mod->str_values[id].fterms.eqrheight, mod->str_values[id].fterms.diamstem, mod->str_values[id].fterms.densstem);
    return temp * c_eda;
  }
  else {
    temp = fr_bedshstr_drag_coef(depth, mod->str_values[id].fterms.eqrheight);
    return temp * c_eda;
  }


}
