function out = model
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');
model.modelPath('C:\Users\Willie D\Desktop\RRAM_model');
%assign parameters
model.param.set('r_f', '3[nm]', 'filament radius');
model.param.set('r_d', '10[nm]', 'device radius');
model.param.set('h_di', '5[nm]', 'height of the dielectric');
model.param.set('h_te', '30[nm]', 'top-electrode height');
model.param.set('h_be', '65[nm]', 'bottom-electrode height');
model.param.set('h_Hf', '10[nm]', 'height Hf');
model.param.set('r_SiO2', '0.5[um]', 'SiO2 radius');
model.param.set('h_SiO2', '0.3[um]', 'SiO2 height');
model.param.set('R_L', '3.1[kohm]', 'Load Resistance');
model.param.set('V_amp', '1.25[V]', 'voltage amplitude');
model.param.set('t_rise', '20[us]', 'pulse rise time');
model.param.set('L_const', '2.44e-8[W/(S*K^2)]');
model.param.set('alpha_fil', '-0.05');
model.param.set('atomic_vibration', '1e-13');
%assign material parameters for TiN/electrode
model.param.set('SHC_TiN', '545.33[J/(kg*K)]', 'Specific Heat Capacity TiN');
model.param.set('Per_TiN', '-1e6', 'Relative Permittivity TiN');
model.param.set('Den_TiN', '5.22e3[kg/m^3]', 'Density TiN');
%assign material parameters for HfO2-x/filament
model.param.set('EC_HfO2x', '7e3[S/m]', 'Electrical Conductivity HfO2-x');
model.param.set('SHC_HfO2x', '140[J/(kg*K)]', 'Specific Heat Capacity HfO2-x');
model.param.set('Per_HfO2x', '-1e6', 'Relative Permittivity HfO2-x');
model.param.set('Den_HfO2x', '12e3[kg/m^3]', 'Density HfO2-x');
%assign material parameter for HfO2/Dielectric
model.param.set('TC_HfO2', '0.5[W/(K*m)]', 'Thermal Conductivity HfO2');
model.param.set('EC_HfO2', '10[S/m]', 'Electrical Conductivity HfO2');
model.param.set('SHC_HfO2', '120[J/(kg*K)]', 'Specific Heat Capacity HfO2');
model.param.set('Per_HfO2', '25', 'Relative Permittivity HfO2');
model.param.set('Den_HfO2', '10e3[kg/m^3]', 'Density HfO2');
%assign material parameter for Hf/metal cap
model.param.set('SHC_Hf', '144[J/(kg*K)]', 'Specific Heat Capacity Hf');
model.param.set('Per_Hf', '-1e6', 'Relative Permittivity Hf');
model.param.set('Den_Hf', '13.31e3[kg/m^3]', 'Density Hf');
%assign material parameter of SiO2
model.param.set('TC_SiO2', '1.38[W/(K*m)]', 'Thermal Conductivity SiO2');
model.param.set('EC_SiO2', '1e-10[S/m]', 'Electrical Conductivity SiO2');
model.param.set('SHC_SiO2', '703[J/(kg*K)]', 'Specific Heat Capacity SiO2');
model.param.set('Per_SiO2', '3.9', 'Relative Permittivity SiO2');
model.param.set('Den_SiO2', '2203[kg/m^3]', 'Density SiO2');

%creates component - local feature
model.component.create('mod1', false);

%use linear extrapolation to use temp dependent TiN electrical resistivity from experiemental data 
model.component('mod1').func.create('int1', 'Interpolation');
model.component('mod1').func('int1').label('TiN resistivity');
model.component('mod1').func('int1').set('table', {'300' '98.74'; '342' '100.85'; '472' '105.76'; '572' '110.73'; '673' '114.95'});
model.component('mod1').func('int1').set('extrap', 'linear');
%use linear extrapolation to use temp dependent Hf electrical resistivity from experiemental data 
model.component('mod1').func.create('int2', 'Interpolation');
model.component('mod1').func('int2').label('Hf');
model.component('mod1').func('int2').set('table', {'293' '33.08';  ...
'300' '34.03';  ...
'350' '40.97';  ...
'400' '48.11';  ...
'450' '55.47';  ...
'500' '63.18';  ...
'550' '70.99';  ...
'600' '78.67';  ...
'650' '86.19';  ...
'700' '93.53';  ...
'750' '100.6';  ...
'800' '107.3';  ...
'850' '113.4';  ...
'900' '119.0';  ...
'950' '124';  ...
'1000' '128.9';  ...
'1100' '138.0';  ...
'1200' '146.4';  ...
'1300' '153.6';  ...
'1400' '159.3';  ...
'1500' '163.5'});
model.component('mod1').func('int2').set('extrap', 'linear');

%creates 2D axisymmetry dimension workspace
model.component('mod1').geom.create('geom1', 2);
model.component('mod1').geom('geom1').axisymmetric(true);
%creates top SiO2
model.component('mod1').geom('geom1').create('r1', 'Rectangle');
model.component('mod1').geom('geom1').feature('r1').label('SiO2 superstrate');
model.component('mod1').geom('geom1').feature('r1').set('pos', {'0' 'h_SiO2 + h_be + h_di + h_Hf +h_te'});
model.component('mod1').geom('geom1').feature('r1').set('size', {'r_SiO2' 'h_SiO2'});
%creates top electrode
model.component('mod1').geom('geom1').create('r2', 'Rectangle');
model.component('mod1').geom('geom1').feature('r2').label('TE');
model.component('mod1').geom('geom1').feature('r2').set('pos', {'0' 'h_SiO2 + h_be + h_di + h_Hf'});
model.component('mod1').geom('geom1').feature('r2').set('size', {'r_d' 'h_te'});
%creates Hf-cap
model.component('mod1').geom('geom1').create('r3', 'Rectangle');
model.component('mod1').geom('geom1').feature('r3').label('Hf-cap');
model.component('mod1').geom('geom1').feature('r3').set('pos', {'0' 'h_SiO2 + h_be + h_di'});
model.component('mod1').geom('geom1').feature('r3').set('size', {'r_d' 'h_Hf'});
%creates filament rectangle
model.component('mod1').geom('geom1').create('r4', 'Rectangle');
model.component('mod1').geom('geom1').feature('r4').label('filament');
model.component('mod1').geom('geom1').feature('r4').set('pos', {'0' 'h_SiO2 + h_be'});
model.component('mod1').geom('geom1').feature('r4').set('size', {'r_f' 'h_di'});
%creates dielectric layer
model.component('mod1').geom('geom1').create('r5', 'Rectangle');
model.component('mod1').geom('geom1').feature('r5').label('Dielectric');
model.component('mod1').geom('geom1').feature('r5').set('pos', {'r_f' 'h_be + h_SiO2'});
model.component('mod1').geom('geom1').feature('r5').set('size', {'r_d-r_f' 'h_di'});
%creates bottom electrode
model.component('mod1').geom('geom1').create('r6', 'Rectangle');
model.component('mod1').geom('geom1').feature('r6').label('BE');
model.component('mod1').geom('geom1').feature('r6').set('pos', {'0' 'h_SiO2'});
model.component('mod1').geom('geom1').feature('r6').set('size', {'r_d' 'h_be'});
%creates SiO2 base
model.component('mod1').geom('geom1').create('r7', 'Rectangle');
model.component('mod1').geom('geom1').feature('r7').label('SiO2 substrate');
model.component('mod1').geom('geom1').feature('r7').set('pos', {'0','0'});
model.component('mod1').geom('geom1').feature('r7').set('size', {'r_SiO2' 'h_SiO2'});
%creates SiO2 enclosure
model.component('mod1').geom('geom1').create('r8', 'Rectangle');
model.component('mod1').geom('geom1').feature('r8').label('SiO2 enclouser');
model.component('mod1').geom('geom1').feature('r8').set('pos', {'r_d' 'h_SiO2'});
model.component('mod1').geom('geom1').feature('r8').set('size', {'r_SiO2-r_d' 'h_be+h_di+h_Hf+h_te'});

%run the geomerty
model.component('mod1').geom('geom1').run;
model.component('mod1').geom('geom1').run('fin');

%define temp and voltage dependendence of filament conductance
model.component('mod1').variable.create('var1');
model.component('mod1').variable('var1').label('filament sigma');
model.component('mod1').variable('var1').set('sigma_fil', 'exp(-alpha_fil*log((t+eps)/(1[s]*atomic_vibration)))*exp(sqrt(abs(V/(8.617e-5[V/K]*T))))');
%model.component('mod1').variable('var1').set('sigma_fil', 'A*exp(sqrt(abs(V/(8.617e-5[V/K]*T))))');
%model.component('mod1').variable('var1').set('sigma_fil', 'exp(-d_E_fil/(8.617e-5[1/K]*T)+sqrt(abs(V/(8.617e-5[V/K]*T))))');
model.component('mod1').variable('var1').selection.geom('geom1', 2);
model.component('mod1').variable('var1').selection.set([3]);

%creates material TiN with the assigned material parameter up above
model.component('mod1').material.create('mat1');
model.component('mod1').material('mat1').label('TiN');
model.component('mod1').material('mat1').propertyGroup('def').set('electricconductivity', '1/(1e-8[ohm*m]*int1(T/1[K]))');
model.component('mod1').material('mat1').propertyGroup('def').set('heatcapacity', 'SHC_TiN');
model.component('mod1').material('mat1').propertyGroup('def').set('relpermittivity','Per_TiN');
model.component('mod1').material('mat1').propertyGroup('def').set('density', 'Den_TiN');
model.component('mod1').material('mat1').propertyGroup('def').set('thermalconductivity', '1/(1e-8[ohm*m]*int1(T/1[K]))*T*L_const');
%creates material HfO2 with the assigned material parameter up above
model.component('mod1').material.create('mat2');
model.component('mod1').material('mat2').label('HfO2');
model.component('mod1').material('mat2').propertyGroup('def').set('electricconductivity','EC_HfO2');
model.component('mod1').material('mat2').propertyGroup('def').set('heatcapacity', 'SHC_HfO2');
model.component('mod1').material('mat2').propertyGroup('def').set('relpermittivity', 'Per_HfO2');
model.component('mod1').material('mat2').propertyGroup('def').set('density', 'Den_HfO2');
model.component('mod1').material('mat2').propertyGroup('def').set('thermalconductivity', 'TC_HfO2');
%creates material HfO2-x with the assigned material parameter up above
model.component('mod1').material.create('mat3');
model.component('mod1').material('mat3').label('HfO2-x');
model.component('mod1').material('mat3').propertyGroup('def').set('electricconductivity', 'EC_HfO2x*sigma_fil');
model.component('mod1').material('mat3').propertyGroup('def').set('heatcapacity', 'SHC_HfO2x');
model.component('mod1').material('mat3').propertyGroup('def').set('relpermittivity', 'Per_HfO2x');
model.component('mod1').material('mat3').propertyGroup('def').set('density', 'Den_HfO2x');
model.component('mod1').material('mat3').propertyGroup('def').set('thermalconductivity', 'EC_HfO2x*sigma_fil*L_const*T');
%creates material SiO2 with the assigned material parameter up above
model.component('mod1').material.create('mat4');
model.component('mod1').material('mat4').label('SiO2');
model.component('mod1').material('mat4').propertyGroup('def').set('electricconductivity', 'EC_SiO2');
model.component('mod1').material('mat4').propertyGroup('def').set('heatcapacity', 'SHC_SiO2');
model.component('mod1').material('mat4').propertyGroup('def').set('relpermittivity', 'Per_SiO2');
model.component('mod1').material('mat4').propertyGroup('def').set('density', 'Den_SiO2');
model.component('mod1').material('mat4').propertyGroup('def').set('thermalconductivity', 'TC_SiO2');
%creates material Hf with the assigned material parameter up above
model.component('mod1').material.create('mat5');
model.component('mod1').material('mat5').label('Hf');
model.component('mod1').material('mat5').propertyGroup('def').set('electricconductivity', '1/(1e-8[ohm*m]*int2(T/1[K]))');
model.component('mod1').material('mat5').propertyGroup('def').set('heatcapacity', 'SHC_Hf');
model.component('mod1').material('mat5').propertyGroup('def').set('relpermittivity', 'Per_Hf');
model.component('mod1').material('mat5').propertyGroup('def').set('density', 'Den_Hf');
model.component('mod1').material('mat5').propertyGroup('def').set('thermalconductivity', '1/(1e-8[ohm*m]*int2(T/1[K]))*T*L_const');

%assign each domain with its corresponding material
model.component('mod1').material('mat1').selection.set([2 5]);
model.component('mod1').material('mat2').selection.set([7]);
model.component('mod1').material('mat3').selection.set([3]);
model.component('mod1').material('mat4').selection.set([1 6 8]);
model.component('mod1').material('mat5').selection.set([4]);

%tells COMOSL to use Electric Current module 
model.component('mod1').physics.create('ec', 'ConductiveMedia', 'geom1');
model.component('mod1').physics('ec').create('gnd1', 'Ground', 1);
model.component('mod1').physics('ec').feature('gnd1').selection.set([4]);
model.component('mod1').physics('ec').create('term1', 'Terminal', 1);
model.component('mod1').physics('ec').feature('term1').selection.set([12]);
model.component('mod1').physics('ec').feature('term1').set('TerminalType', 'Circuit');
%model.component('mod1').physics('ec').create('ci1', 'ContactImpedance', 1);
%model.component('mod1').physics('ec').feature('ci1').selection.set([6 8 15 16]);
%model.component('mod1').physics('ec').feature('ci1').set('spec_type', 'surfimp');
%model.component('mod1').physics('ec').feature('ci1').set('rhos', 'rho_cont*exp(-rho_act/(8.617e-5[1/K]*T))');
%tells COMSOL to use heat transfer module
model.component('mod1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('mod1').physics('ht').prop('AmbientSettings').set('T_amb', '298.15[K]');
model.component('mod1').physics('ht').feature('init1').set('Tinit_src', 'root.mod1.ht.T_amb');
model.component('mod1').physics('ht').create('temp1', 'TemperatureBoundary', 1);
model.component('mod1').physics('ht').feature('temp1').set('T0_src', 'root.mod1.ht.T_amb');
model.component('mod1').physics('ht').feature('temp1').selection.set([2 13]);
%heat can escape in radiative process
model.component('mod1').physics('ht').create('ds1', 'DiffuseSurface', 1);
model.component('mod1').physics('ht').feature('ds1').set('Tamb_src', 'root.mod1.ht.T_amb');
model.component('mod1').physics('ht').feature('ds1').selection.set([4 6 8 10 12 14 15 16 17 18 19 20 21 22]);
model.component('mod1').physics('ht').feature('ds1').set('epsilon_rad', 0.9);
model.component('mod1').physics('ht').feature('ds1').set('epsilon_rad_mat', 'userdef');

%tells COMSOL to use electrical circuit module
model.component('mod1').physics.create('cir', 'Circuit', 'geom1');
%set ramping voltage pulse as voltage source
model.component('mod1').physics('cir').create('V1', 'VoltageSource', -1);
model.component('mod1').physics('cir').feature('V1').set('Connections', [1; 0]);
model.component('mod1').physics('cir').feature('V1').set('sourceType', 'PulseSource');
model.component('mod1').physics('cir').feature('V1').set('value', 'V_amp');
model.component('mod1').physics('cir').feature('V1').set('td', 't_rise/2');
model.component('mod1').physics('cir').feature('V1').set('tr', 't_rise');
model.component('mod1').physics('cir').feature('V1').set('tf', 't_rise');
model.component('mod1').physics('cir').feature('V1').set('pw', 't_rise');
model.component('mod1').physics('cir').feature('V1').set('Tper', '4*t_rise');
%connect load resistor in series
model.component('mod1').physics('cir').create('R1', 'Resistor', -1);
model.component('mod1').physics('cir').feature('R1').set('R', 'R_L');
model.component('mod1').physics('cir').feature('R1').set('Connections', [1; 2]);
%connect the circuit with device on terminal defined above
model.component('mod1').physics('cir').create('IvsU1', 'ModelDeviceIV', -1);
model.component('mod1').physics('cir').feature('IvsU1').set('Connections', [2; 0]);
model.component('mod1').physics('cir').feature('IvsU1').set('V_src', 'root.mod1.ec.V0_1');
%couples electric current model and heat transfer module
model.component('mod1').multiphysics.create('emh1', 'ElectromagneticHeating', -1);
model.component('mod1').multiphysics.create('tc1', 'TemperatureCoupling', -1);

%creates mesh
model.component('mod1').mesh.create('mesh1');
model.component('mod1').mesh('mesh1').feature('size').set('hauto', 1);
%fine mesh for SiO2
model.component('mod1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('mod1').mesh('mesh1').feature('ftri1').label('SiO2');
model.component('mod1').mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.component('mod1').mesh('mesh1').feature('ftri1').selection.set([1 6 8]);
model.component('mod1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', 4);
%extra fine mesh for electrodes+Hf-cap
model.component('mod1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('mod1').mesh('mesh1').feature('ftri2').label('electrode');
model.component('mod1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('mod1').mesh('mesh1').feature('ftri2').selection.set([2 4 5]);
model.component('mod1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('ftri2').feature('size1').set('hauto', 2);
%custom mesh for insulator layer - max mesh size allowed - 1.5nm
model.component('mod1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('mod1').mesh('mesh1').feature('ftri3').label('insulator');
model.component('mod1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('mod1').mesh('mesh1').feature('ftri3').selection.set([3 7]);
model.component('mod1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('ftri3').feature('size1').set('hauto', 1);
model.component('mod1').mesh('mesh1').feature('ftri3').feature('size1').set('custom', 'on');
model.component('mod1').mesh('mesh1').feature('ftri3').feature('size1').set('hmax', 'r_f');
model.component('mod1').mesh('mesh1').feature('ftri3').feature('size1').set('hmaxactive', true);
model.component('mod1').mesh('mesh1').run;

%creates a time-dependent transient solver and solves the above problem
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tunit', 'ms');
model.study('std1').feature('time').set('tlist', 'range(0,t_rise/30,3/2*t_rise)');
model.study('std1').run;
out = model;
