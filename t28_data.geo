
//-----your variables-----

mm=1e-3;
//height of the coils is 52 mm


//---important variables----

width_al = 1*mm;//width of the AirLayer in [m]
mass = 0.107;// mass of the plate in [kg]
gravity_const = 9.8;
Rri = 140*mm;
Rro = 160*mm;
lr0 = 5*mm;
lal = 1*mm;

//Coil parameters
width_ic = 28*mm;
width_oc = 15*mm;
height_coil = 52*mm;
ic_radius = 41*mm;
oc_radius = ic_radius + 46.5*mm;
lci = 1*mm;
lco = 1*mm;

//Plate parameters
plate_radius = 65*mm;
plate_thick = 3*mm;
initial_height = 3.8*mm;
lcp = 0.5*mm;

//------display variables----------------

ppm = "Input/11Mesh control (Nbr of divisions)/";
colorpp = "Ivory";
colorro  = "LightGrey";
close_menu = 0;
DefineConstant[
  md = { 1., Name StrCat[ppm, "0Mesh density"],
    Highlight Str[colorpp], Closed close_menu},
  nn_hplate   = { Ceil[md*2], Name StrCat[ppm, "10Plate width"], ReadOnly 1,
    Highlight Str[colorro], Closed close_menu},
  nn_coil1  = { Ceil[md*5], Name StrCat[ppm, "11Inner coil width"], ReadOnly 1,
    Highlight Str[colorro]},
  nn_coil2  = { Ceil[md*3], Name StrCat[ppm, "11Inner coil width"], ReadOnly 1,
    Highlight Str[colorro]},
  nn_ri = { Ceil[md*8], Name StrCat[ppm, "2One fourth shell in"], ReadOnly 1,
    Highlight Str[colorro]},
  nn_ro = { Ceil[md*6], Name StrCat[ppm, "3One fourth shell out"], ReadOnly 1,
    Highlight Str[colorro]}
];


color_time="AliceBlue";
color_time_var="Ivory";
color_mec ="Ivory";
color_ele ="Ivory";


input_mec="Input/Mechanics/";
output_mec="Output/Mechanics/";
output_mag="Output/Magnetics/";

output="Output/";

str_vel = "3v_y [m\s]"; //3Vertical velocity [m\s]
str_disp = "2d_y [m]";  // 2Vertical displacement [m]
ExtGnuplot  = ".dat";
ExtGmsh     = ".pos";

M_calc = StrCat[output,"01to post-process/"] ;

DefineConstant[
  Flag_AnalysisType = {1,  Choices{0="Frequency domain",  1="Time domain"},
    Name "Input/00Type of analysis", Highlight "Blue",
    Help Str["- Use 'Frequency domain' to compute the steady-state fields at equilibrium position",
      "- Use 'Time domain' to compute the dynamic response of the conductive plate"]}
];

p_init = (Flag_AnalysisType==0) ? 11.3*1e-3 : 3.8*1e-3 ; // plate initial position

If(Flag_AnalysisType==0)// Steady state - equilibrium at 11.3 mm (measurements)
  DefineConstant[ displacementY = { p_init,
      Name StrCat["Input/",str_disp], Visible 0} ];
  UndefineConstant[ "Input/1step" ];
  UndefineConstant[ StrCat[output_mec, str_disp] ];
EndIf
If(Flag_AnalysisType==1) // transient
  DefineConstant[ displacementY = { p_init, Min -p_init, Max 100*1e-3,
      Name StrCat[output_mec, str_disp], ReadOnlyRange 1, Highlight Str[color_mec]} ];
  UndefineConstant[ StrCat["Input/",str_disp] ];
EndIf

INNER_COIL = 1000;
OUTER_COIL = 2000;
PLATE = 3000;
AIR_SHELL = 4000;
INNER_SHELL = 5000;
AIR_LAY_UP = 6000;
AIR_LAY_RIGHT = 7000;
AIR_LAY_DOWN = 8000;
OUTER_BND = 9000;
AXIS_SYM = 10000;
SKIN_DOMAIN = 11000;
