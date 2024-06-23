mlexp (choice_x*ln(exp({a}*((1-{b}*s_x-{c}*r_x-{d}*q-{e}*v)*self_x^{f}+({b}*s_x+{c}*r_x+{d}*q+{e}*v)*other_x^{f})^(1/{f})) /*
*/ /(exp({a}*((1-{b}*s_x-{c}*r_x-{d}*q-{e}*v)*self_x^{f}+({b}*s_x+{c}*r_x+{d}*q+{e}*v)*other_x^{f})^(1/{f}))+exp({a}*((1-{b}*s_y-{c}*r_y-{d}*q-{e}*v)*self_y^{f}+({b}*s_y+{c}*r_y+{d}*q+{e}*v)*other_y^{f})^(1/{f})))) /*
*/ +(1-choice_x)*ln(exp({a}*((1-{b}*s_y-{c}*r_y-{d}*q-{e}*v)*self_y^{f}+({b}*s_y+{c}*r_y+{d}*q+{e}*v)*other_y^{f})^(1/{f})) /*
*/ /(exp({a}*((1-{b}*s_x-{c}*r_x-{d}*q-{e}*v)*self_x^{f}+({b}*s_x+{c}*r_x+{d}*q+{e}*v)*other_x^{f})^(1/{f}))+exp({a}*((1-{b}*s_y-{c}*r_y-{d}*q-{e}*v)*self_y^{f}+({b}*s_y+{c}*r_y+{d}*q+{e}*v)*other_y^{f})^(1/{f}))))) /*
*/, vce(cluster sid) from(a=0.13 b=0.06 c=0.25 d=0.07 e=-0.05 f=1)



mlexp (choice_a*ln(exp({a}*(1-{b}*lambdas_a))))



import delimited "C:\Users\User\OneDrive\Documents\Academics\Data\Bruhin2019\choices_exp1_inclZalloc_inclSubjChar.csv", clear
drop self_z other_z
rename s_x ss_x
rename s_y ss_y
rename r_y rr_y
rename r_x rr_x
rename q qq
rename v vv
gen id = _n
gen choice_y = 1 - choice_x
reshape long self other choice ss rr, i(id) j(alt) string

mlexp (choice*ln(exp({a}*((1-{b}*ss_x-{c}*r_x-{d}*q-{e}*v)*self_x^{f}+({b}*s_x+{c}*r_x+{d}*q+{e}*v)*other_x^{f})^(1/{f})) /*
*/ /(exp({a}*((1-{b}*s_x-{c}*r_x-{d}*q-{e}*v)*self_x^{f}+({b}*s_x+{c}*r_x+{d}*q+{e}*v)*other_x^{f})^(1/{f}))+exp({a}*((1-{b}*s_y-{c}*r_y-{d}*q-{e}*v)*self_y^{f}+({b}*s_y+{c}*r_y+{d}*q+{e}*v)*other_y^{f})^(1/{f})))) /*