# objects=param_mod.o leopard_rg.o disp_det_rg.o Z_func.o dZ_func.o integrator.o muller_rg.o polyfit.o read_data_rg.o read_distr_rg.o spline_interpol.o exp_Bessel_In.o exp_dBessel_In.o cerror.o cont_frac.o Bessel_int.o int_para.o int_para_mpfun.o F12.o F12_mpfun.o F23.o F23_mpfun.o gamma_func.o acc_F.o acc_Kvpa.o get_splinecoeff.o exp_Bessel_In_mpfun.o fort_Bes.o
objects=param_mod.o leopardv3_mod.o leopard_rg.o disp_det_rg.o Z_func.o dZ_func.o integrator.o polyfit.o read_data_rg.o read_distr_rg.o spline_interpol.o exp_Bessel_In.o exp_dBessel_In.o cerror.o cont_frac.o Bessel_int.o int_para.o int_para_mpfun.o F12.o F12_mpfun.o F23.o F23_mpfun.o gamma_func.o acc_F.o acc_Kvpa.o get_splinecoeff.o exp_Bessel_In_mpfun.o fort_Bes.o
mpfun_obj = mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq2.o  mpmodule.o
mpfun_mod = mpfuna.mod  mpfunb.mod mpfunc.mod mpfund.mod mpfune.mod mpfunf.mod mpfung.mod  mpmodule.mod
f90comp = gfortran
options =  -fdefault-real-8 -O3 -g -ffpe-trap=invalid -ffpe-trap=zero -ffpe-trap=overflow
options_mp = -O3 -g

leopard_rg: $(objects) $(mpfun_obj)
	$(f90comp) -o leopard_rg $(options) $(objects) $(mpfun_obj)

param_mod.mod: param_mod.o param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

param_mod.o: param_mod.f90
	$(f90comp) -c $(options) param_mod.f90

leopardv3_mod.mod: leopardv3_mod.o leopardv3_mod.f90
	$(f90comp) -c $(options) leopardv3_mod.f90

leopardv3_mod.o: leopardv3_mod.f90
	$(f90comp) -c $(options) leopardv3_mod.f90

mpmodule.mod: mpmodule.o mpmodule.f90
	$(f90comp) -c $(options_mp) mpmodule.f90

mpmodule.o: mpmodule.f90
	$(f90comp) -c $(options_mp) mpmodule.f90

mpfuna.mod: mpfuna.o mpfuna.f90
	$(f90comp) -c $(options_mp) mpfuna.f90

mpfuna.o: mpfuna.f90
	$(f90comp) -c $(options_mp) mpfuna.f90

mpfunb.mod: mpfunbq.o mpfunbq.f90
	$(f90comp) -c $(options_mp) mpfunbq.f90

mpfunbq.o: mpfunbq.f90
	$(f90comp) -c $(options_mp) mpfunbq.f90

mpfunc.mod: mpfunc.o mpfunc.f90
	$(f90comp) -c $(options_mp) mpfunc.f90

mpfunc.o: mpfunc.f90
	$(f90comp) -c $(options_mp) mpfunc.f90

mpfund.mod: mpfund.o mpfund.f90
	$(f90comp) -c $(options_mp) mpfund.f90

mpfund.o: mpfund.f90
	$(f90comp) -c $(options_mp) mpfund.f90

mpfune.mod: mpfune.o mpfune.f90
	$(f90comp) -c $(options_mp) mpfune.f90

mpfune.o: mpfune.f90
	$(f90comp) -c $(options_mp) mpfune.f90

mpfunf.mod: mpfunf.o mpfunf.f90
	$(f90comp) -c $(options_mp) mpfunf.f90

mpfunf.o: mpfunf.f90
	$(f90comp) -c $(options_mp) mpfunf.f90

mpfung.mod: mpfungq2.o mpfungq2.f90
	$(f90comp) -c $(options_mp) mpfungq2.f90

mpfungq2.o: mpfungq2.f90
	$(f90comp) -c $(options_mp) mpfungq2.f90

leopard_rg.o: param_mod.mod leopard_rg.f90
	$(f90comp) -c $(options) leopard_rg.f90

disp_det_rg.o: param_mod.mod disp_det_rg.f90
	$(f90comp) -c $(options) disp_det_rg.f90

integrator.o: param_mod.mod integrator.f90
	$(f90comp) -c $(options)  integrator.f90

Z_func.o: param_mod.mod Z_func.f90
	$(f90comp) -c $(options)  Z_func.f90

dZ_func.o: dZ_func.f90
	$(f90comp) -c $(options)  dZ_func.f90

cerror.o: param_mod.mod cerror.f90
	$(f90comp) -c $(options)  cerror.f90

cont_frac.o: cont_frac.f90
	$(f90comp) -c $(options)  cont_frac.f90

# muller_rg.o: param_mod.mod muller_rg.f90
# 	$(f90comp) -c $(options)  muller_rg.f90

polyfit.o: polyfit.f90
	$(f90comp) -c $(options)  polyfit.f90

read_data_rg.o: param_mod.mod read_data_rg.f90
	$(f90comp) -c $(options)  read_data_rg.f90

exp_Bessel_In.o: param_mod.mod  exp_Bessel_In.f90
	$(f90comp) -c $(options)  exp_Bessel_In.f90

exp_Bessel_In_mpfun.o: param_mod.mod  $(mpfun_mod) exp_Bessel_In_mpfun.f90
	$(f90comp) -c $(options)  exp_Bessel_In_mpfun.f90

exp_dBessel_In.o: exp_dBessel_In.f90
	$(f90comp) -c $(options)  exp_dBessel_In.f90

Bessel_int.o: param_mod.mod Bessel_int.f90
	$(f90comp) -c $(options)  Bessel_int.f90

read_distr_rg.o: param_mod.mod read_distr_rg.f90
	$(f90comp) -c $(options)  read_distr_rg.f90

spline_interpol.o: spline_interpol.f90
	$(f90comp) -c $(options)  spline_interpol.f90

get_splinecoeff.o: param_mod.mod get_splinecoeff.f90
	$(f90comp) -c $(options)  get_splinecoeff.f90

F23.o: F23.f90
	$(f90comp) -c $(options)  F23.f90

F23_mpfun.o: $(mpfun_mod) F23_mpfun.f90
	$(f90comp) -c $(options)  F23_mpfun.f90

F12.o: F12.f90
	$(f90comp) -c $(options)  F12.f90

F12_mpfun.o: $(mpfun_mod) F12_mpfun.f90
	$(f90comp) -c $(options)  F12_mpfun.f90

acc_F.o: acc_F.f90
	$(f90comp) -c $(options)  acc_F.f90

acc_Kvpa.o: acc_Kvpa.f90
	$(f90comp) -c $(options)  acc_Kvpa.f90

int_para.o: param_mod.mod int_para.f90
	$(f90comp) -c $(options)  int_para.f90

int_para_mpfun.o: param_mod.mod $(mpfun_mod) int_para_mpfun.f90
	$(f90comp) -c $(options)  int_para_mpfun.f90

gamma_func.o: gamma_func.f90
	$(f90comp) -c $(options)  gamma_func.f90

fort_Bes.o: fort_Bes.f90
	$(f90comp) -c $(options)  fort_Bes.f90

clean:
	rm leopard_rg 
	rm param_mod.mod
	rm leopardv3_mod.mod
	rm $(objects)
	rm $(mpfun_mod)
	rm $(mpfun_obj)
