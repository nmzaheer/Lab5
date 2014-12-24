Include "t28_data.geo";

ResDir="resultaten"; // You may give here a name of a directory, e.g. "res/", results will then be saved there

DefineConstant[
  Clean_Results = {1, Name "Input/01Remove previous result files", Choices {0,1}},
  time_min = { 0, Visible (Flag_AnalysisType == 1),
               Name "Input/30time min.", Highlight Str[color_time]},
  time_max = { 1000e-3, Min time_min, Visible (Flag_AnalysisType == 1),
               Name "Input/31time max.", Highlight Str[color_time]},
  delta_time = { 0.2e-3, Min time_max/5000, Visible (Flag_AnalysisType == 1),
                 Name "Input/32delta time", Highlight Str[color_time]},
  NbSteps = (time_max - time_min)/delta_time,


  // Do not change the step constant ==> it controls the time loop inside onelab
   step = {0, Min 0, Max NbSteps, Step 1, Loop  (Flag_AnalysisType == 1),
          Name "Input/34step", Visible (Flag_AnalysisType == 1), Highlight Str[color_time_var]},
  Flag_Cir = 0,
  Flag_NL = { 0, Choices{0,1}, Name "Input/60Nonlinear BH-curve", Visible 0},
  
  CoefGeo = 2*Pi //Axisymmetry
];	


Group {
  // Define the Physical groups taking into account your geometry:
  // List of groups appearing in the formulation that you should complete
    
    DomainB= Region[{INNER_COIL,OUTER_COIL}]; //domain with all stranded coils (allowing for circuit coupling), 'B' of 'bobiné'
    DomainC_Moving= Region[ {PLATE}];//domain of the conducting (~metals) material that can move
    DomainC_NonMoving=Region [{}];//domain of the conducting (~metals) material that can not move. no coils
    DomainC = Region [{DomainC_Moving, DomainC_NonMoving}];//domain of all conducting materials
    DomainCC= Region [{AIR_SHELL,INNER_SHELL,AIR_LAY_UP,AIR_LAY_RIGHT,AIR_LAY_DOWN, DomainB}];//domain of 'complementary conducting' materials, everything that is not in Domain C
    DomainInf=Region [{}]; // Only used when applying a transformation to infinity
	
	SkinDomainC_Moving=Region[{SKIN_DOMAIN}]; // The outer skin of the DomainC_Moving of moving conducting materials
	AirLayer=Region[{AIR_LAY_UP,AIR_LAY_RIGHT,AIR_LAY_DOWN}]; //Airlayer around the moving material
	AirTot = Region[{AIR_SHELL,INNER_SHELL,AIR_LAY_UP,AIR_LAY_RIGHT,AIR_LAY_DOWN}]; //Collection of all Air elements
	SurfaceGInf = Region[{OUTER_BND}]; //Boundary at infinity
	DomainS = Region[{}];//empty
	
	
	Coil_Inside=Region[{INNER_COIL}];//collection of the inner coil(s) 
	Coil_Outside=Region[{OUTER_COIL}];//collection of the outer coil(s)
	Coils = Region[{INNER_COIL,OUTER_COIL}]; //collection of all coils
		
	SymmetryAxis = Region[{AXIS_SYM}];//Axis of symmetry 
  
  // Do not remove these two domains
  DomainKin = #1234 ; // Dummy region number for mechanical equation
  DomainDummy = #12345 ; // Dummy region number for postpro with functions
  
  DomainNL = Region[{}]; //Domain of Non-Linear Material
  DomainL  = Region[{ DomainC, DomainCC }]; //Domain of Linear Material
  Domain = Region[ {DomainL, DomainNL} ] ; // Domain of all Material
  
  
  DefineGroup[Dummy, E1, E2, DomainZt_Cir]; //not used
}



Function {
  DefineFunction[ js, dhdb_NL ];

  mu0 = 4e-7*Pi ;
  mur = 1.0 ;
  nu[#{Coils, AirTot, DomainC_Moving}]  = 1. / mu0 ; //definition of nu

  // change for non-linear calculations
  Nb_max_iter        = 30 ; // Number of Maximal iterations. default value=30;
  relaxation_factor  = 1; //relaxation factor, default value=1
  stop_criterion     = 1e-5; //min error at which to stop. default value=1e-5

  //Val_Rint = ***; //only used when applying infinity transformation
  //Val_Rext = ***;

  sigma[#{DomainC_Moving}] = 3.4e7; //aluminium
  sigma[#{Coils}] = 5.96e7; //cupper
  i0 = 20;

  NbWires[#{Coil_Inside}] = 960;
  NbWires[#{Coil_Outside}] = 576;
  SurfCoil[#{Coil_Inside}] = SurfaceArea[]{INNER_COIL}; //surface of inside coil(s). In [m]
  SurfCoil[#{Coil_Outside}] = SurfaceArea[]{OUTER_COIL}; //surface of outside coil(s). In [m]
  Idir[#{Coil_Inside}] =  1.0; // Direction of current in inside coil. Scalar value, e.g= -1.0
  Idir[#{Coil_Outside}] = 1.0; // Direction of current in outside coil. Scalar value, e.g= -1.0

  DefineConstant[
    Freq = { 50, Name "Input/50Frequency [Hz]", Highlight "AliceBlue"},
    II = { i0,  Name "Input/51Current [A]", Highlight "AliceBlue"},
    velocityY = { 0., Name StrCat[output_mec, str_vel], // Vertical velocity
      Visible (Flag_AnalysisType == 1),  Highlight Str[color_mec]}
    forceY = { 0., Name StrCat[output_mag,"41F_y MST [N]"], // Attraction force, for updating graph correctly
      Visible 1,  Highlight Str[color_ele]}
  ];

  I1[] = F_Sin_wt_p[]{2*Pi*Freq, (Flag_AnalysisType==1)?0:Pi/2} ;
  I2[] = F_Sin_wt_p[]{2*Pi*Freq, (Flag_AnalysisType==1)?0:Pi/2} ;

  
  
// HINT: to find the filling factor of the windings, you need the number of wires and the crossection of 1 wire 
// with respect to the total crossection of the coil.
// you can find the crossection of 1 wire by calculating back from the American Wire Gauge-standard. 18 AWG ==> diameter ==> crossection
  FillFactor_Winding[#{Coil_Inside}] = 0.5407;
  FillFactor_Winding[#{Coil_Outside}] = 0.6055;

  Rb[] = CoefGeo*FillFactor_Winding[]*NbWires[]^2/SurfCoil[]/sigma[]; //look in course notes 
  Resistance[] = Rb[];

  time0 = time_min + step * delta_time;

  dXYZ[#{DomainC_NonMoving, DomainB}] = Vector[0., 0., 0.];
  dXYZ[#{DomainC_Moving}] = Vector[0, DefineNumber[0, Name "DeltaU", Visible 0], 0];
  a_previousstep[] = Vector[0, 0, Field[XYZ[]-dXYZ[]]] ;

  // Normal for computing Force with Maxwell stress tensor
  p_current = displacementY ; // Current position
  N[#{AirLayer}] = 1/width_al *
  ( (Y[] >= 3*1e-3+p_current) ? Vector[ 0.,  1.,0.] :
    (Y[] <= p_current) ? Vector[ 0., -1.,0.] :
    (X[] >= 130*1e-3/2)  ? Vector[ 1.,  0.,0.] : Vector[ -1., 0.,0.] ) ;

  // Maxwell tensor
  Tmax[] = (SquDyadicProduct[$1]-SquNorm[$1]*TensorDiag[0.5,0.5,0.5])/mu0 ;
  Tmax_cplx[] = Re[0.5*(TensorV[CompX[$1]*Conj[$1], CompY[$1]*Conj[$1], CompZ[$1]*Conj[$1]] - $1*Conj[$1]*TensorDiag[0.5, 0.5, 0.5])/mu0] ;
  Tmax_cplx_2f[] = 0.5*(TensorV[CompX[$1]*$1, CompY[$1]*$1, CompZ[$1]*$1] - $1*$1 * TensorDiag[0.5, 0.5, 0.5])/mu0 ;// TO CHECK...

  // Kinematics -> Movement
  Inertia[] = mass ;
  Fmg[] = mass * gravity_const;
  // Using Registers for storing values, previous post-processing is needed
  Fmag[] = CompY[#55] ;
  Fmag_vw[] = CompY[#56] ;
  Florentz[] = CompY[#57] ;

  Fmag_FD[] = CompY[#65] ;
  Fmag_vw_FD[] = CompY[#67] ;
  Florentz_FD[] = CompY[#66] ;
}

Jacobian {
  // Axisymmetry
  { Name Vol ;
    Case { //{ Region DomainInf ;
      //  Jacobian VolAxiSquSphShell {Val_Rint, Val_Rext} ;}
      { Region All ; Jacobian VolAxiSqu ; }
    }
  }
  { Name Sur ;
    Case { { Region All ; Jacobian SurAxi ; }
    }
  }
}

Integration {
  { Name I1 ; Case { { Type Gauss ; Case {
          { GeoElement Triangle   ; NumberOfPoints  6 ; }
          { GeoElement Quadrangle ; NumberOfPoints  4 ; }
          { GeoElement Line       ; NumberOfPoints  13 ; }
        } }
    }
  }
}


Constraint {
  { Name MVP_2D ;
    Case {
      { Region SurfaceGInf ; Type Assign; Value 0. ; }
      { Region SymmetryAxis ; Type Assign; Value 0. ; }

      If(time0 != 0 && Flag_AnalysisType==1)
        { Region Domain ; Type InitFromResolution ; NameOfResolution ProjectionInit ; }
      EndIf
    }
  }

  { Name Current_2D ;
    Case {
      { Region Coil_Inside ; Value II*Idir[] ; TimeFunction  I1[] ;}
      { Region Coil_Outside ; Value II*Idir[] ; TimeFunction  I2[] ;}
    }
  }
  { Name Voltage_2D ;
    Case {
    }
  }

  { Name Current_Cir ;
    Case {
    }
  }
  { Name Voltage_Cir ;
    Case {
      //{ Region E1 ; Value VV1 ; TimeFunction  Voltage1[] ;}
      //{ Region E2 ; Value VV2 ; TimeFunction  Voltage2[] ;}
    }
  }

  { Name ElectricalCircuit ; Type Network ;
    Case Circuit1 {
     { Region E1    ; Branch {100,102} ; }
     { Region Coil_Inside ; Branch {102,100} ; }
    }
    Case Circuit2 {
     { Region E2    ; Branch {200,202} ; }
     { Region Coil_Outside ; Branch {202,200} ; }
    }
  }

  { Name CurrentPosition ;
    Case {
      { Region DomainKin ; Type Init ; Value displacementY; }
    }
  }
  { Name CurrentVelocity ;
    Case {
      { Region DomainKin ; Type Init ; Value velocityY ; }
    }
  }

}

FunctionSpace {
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se1 ; NameOfCoef ae1 ; Function BF_PerpendicularEdge ;
        Support Domain ; Entity NodesOf [ All ] ; }
    }
    Constraint {
      { NameOfCoef ae1 ; EntityType NodesOf  ; NameOfConstraint MVP_2D ; }//Error   : GetDP - Non vector argument for function 'CompY'

    }
  }

    { Name Hregion_i_Mag_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainB ; Entity DomainB ; }
    }
    GlobalQuantity {
      { Name Ib ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Ub ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Ub ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Ib ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }

  { Name Hregion_Z ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainZt_Cir ; Entity DomainZt_Cir ; }
    }
    GlobalQuantity {
      { Name Iz ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Uz ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Uz ; EntityType Region ; NameOfConstraint Voltage_Cir ; }
      { NameOfCoef Iz ; EntityType Region ; NameOfConstraint Current_Cir ; }
    }
  }

  { Name Position ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainKin ; Entity DomainKin ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf  ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType Region ; NameOfConstraint CurrentPosition ; }
    }
  }

  { Name Velocity ; Type Scalar ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_Region ;
        Support DomainKin ; Entity DomainKin ; } }
    GlobalQuantity {
      { Name V ; Type AliasOf  ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef V ; EntityType Region ; NameOfConstraint CurrentVelocity ; }
    }
  }

}

Formulation {

 { Name Projection ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
    }
    Equation {
      Galerkin { [  Dof{a}, {a} ] ;
        In Region[{DomainC, DomainB}] ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ -a_previousstep[], {a} ] ;
        In Region[{DomainC, DomainB}] ; Jacobian Vol ; Integration I1 ; }
    }
  }

  { Name MagDyn_a_2D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
      { Name ir ; Type Local  ; NameOfSpace Hregion_i_Mag_2D ; }

      { Name Ub ; Type Global ; NameOfSpace Hregion_i_Mag_2D [Ub] ; }
      { Name Ib ; Type Global ; NameOfSpace Hregion_i_Mag_2D [Ib] ; }
      { Name Uz ; Type Global ; NameOfSpace Hregion_Z [Uz] ; }
      { Name Iz ; Type Global ; NameOfSpace Hregion_Z [Iz] ; }
    }
    Equation {
      Galerkin { [ nu[{d a}] * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }
      Galerkin { JacNL[ dhdb_NL[{d a}] * Dof{d a} , {d a} ]  ;
        In DomainNL ; Jacobian Vol ; Integration I1 ; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -js[], {a} ] ;
        In DomainS ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires[]/SurfCoil[] * Dof{ir} , {a} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ NbWires[]/SurfCoil[] * Dof{a} , {ir} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }

      GlobalTerm { [ Resistance[]  * Dof{Ib}, {Ib} ] ; In DomainB ;}
      GlobalTerm { [ Dof{Ub}/CoefGeo , {Ib} ] ; In DomainB ; }

      If(Flag_Cir)
        GlobalTerm { [ Dof{Uz}        , {Iz} ] ; In Resistance_Cir ; }
        GlobalTerm { [ Resistance[{Iz}] * Dof{Iz} , {Iz} ] ; In Resistance_Cir ; }

        GlobalEquation {
          Type Network ; NameOfConstraint ElectricalCircuit ;
          { Node {Ib}; Loop {Ub}; Equation {Ub}; In DomainB ; }
          { Node {Iz}; Loop {Uz}; Equation {Uz}; In DomainZt_Cir ; }
        }
      EndIf
    }
  }

  { Name Mechanical ; Type FemEquation ;
    Quantity {
      { Name V ; Type Global ; NameOfSpace Velocity [V] ; } // velocity
      { Name U ; Type Global ; NameOfSpace Position [U] ; } // position
    }
    Equation {
      GlobalTerm { DtDof [ Inertia[] * Dof{V} , {V} ] ; In DomainKin ; }
      GlobalTerm { NeverDt[   Fmg[] , {V} ] ; In DomainKin ; }
      GlobalTerm { NeverDt[ -Fmag[] , {V} ] ; In DomainKin ; }

      GlobalTerm { DtDof [ Dof{U} , {U} ] ; In DomainKin ; }
      GlobalTerm {       [-Dof{V} , {U} ] ; In DomainKin ; }
    }
  }

}

Resolution {
  { Name ProjectionInit ;
    System {
      { Name Pr ; NameOfFormulation Projection ; DestinationSystem A; }
    }
    Operation {
      GmshRead[ StrCat[ResDir,"tmp.pos"] ]; Generate[Pr] ; Solve[Pr] ; TransferSolution[Pr] ;
    }
  }

  { Name  Analysis ;
    System {
      If(Flag_AnalysisType ==0)
        { Name A ; NameOfFormulation MagDyn_a_2D ; Type Complex; Frequency Freq; }
      EndIf
      If(Flag_AnalysisType ==1)
        { Name A ; NameOfFormulation MagDyn_a_2D ; }
        { Name M ; NameOfFormulation Mechanical ; }
      EndIf
    }
    Operation {
      CreateDir[ResDir];
      If(Clean_Results==1 && (step==0))
        DeleteFile[
          StrCat[ResDir, "Fmag_*", ExtGnuplot]
        ];
        DeleteFile[StrCat[ResDir, "Fmag_vw", ExtGnuplot]];
        DeleteFile[StrCat[ResDir, "Fmag_vw_cplx", ExtGnuplot]];
        DeleteFile[StrCat[ResDir, "Florentz", ExtGnuplot]];
        DeleteFile[StrCat[ResDir, "disp", ExtGnuplot]];
        DeleteFile[StrCat[ResDir, "velocity", ExtGnuplot]];
      EndIf
      If(Flag_AnalysisType ==0)
          Generate[A]; Solve[A];
          SaveSolution[A];
          PostOperation[MapMag_freq] ;
      EndIf

      If(Flag_AnalysisType ==1)
        SetTime[time0] ;
        InitSolution[A] ; SaveSolution[A] ;
        // Init M twice for avoiding warning (a = d_t^2 x)
        InitSolution[M] ; InitSolution[M]; SaveSolution[M] ;

        TimeLoopTheta[time0, time0+delta_time, delta_time, 1.]{ // do 1 step
          If(!Flag_NL)
            Generate[A]; Solve[A];
          EndIf
          If(Flag_NL)
            IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor] {
              GenerateJac[A] ; SolveJac[A] ; }
          EndIf
          SaveSolution[A];
          PostOperation[MapMag_min] ;

          Generate[M] ; Solve[M] ; SaveSolution[M] ;
          PostOperation[MapMec] ;
        }
      EndIf
    }
  }
EndIf

}


PostProcessing {

 { Name MagDyn_a_2D ; NameOfFormulation MagDyn_a_2D ;
   PostQuantity {
     { Name a ; Value { Term { [  {a}*X[] ]   ; In Domain ; Jacobian Vol ; } } }
     { Name az ; Value { Term { [  CompZ[{a}]*X[] ]   ; In Domain ; Jacobian Vol ; } } } //Axisymmetry
     { Name az_nox ; Value { Term { [  CompZ[{a}] ]   ; In Domain ; Jacobian Vol ; } } } //Needed for proyection

     { Name b  ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
     { Name j  ; Value { Term { [ -sigma[] * Dt[{a}] ] ; In DomainC ; Jacobian Vol ; } } }
     { Name jz  ; Value { Term { [ CompZ[-sigma[] * Dt[{a}]] ] ; In Domain ; Jacobian Vol ; } } }

     { Name boundary  ; Value { Term { [ 1 ] ; In Dummy ; Jacobian Vol ; } } } // Dummy quantity - for visualization

     { Name Flux ; Value { Integral { [ CoefGeo *Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
           In Coils  ; Jacobian Vol ; Integration I1 ; } } }
     { Name W  ; Value { Integral { [ CoefGeo *nu[]/2*SquNorm[{d a}] ] ;
           In Domain ; Jacobian Vol ; Integration I1 ; } } }

     { Name F  ; Value { Integral { [ CoefGeo * Tmax[{d a}] * N[] ] ;
           In Domain ; Jacobian Vol ; Integration I1 ; } } }
	 { Name Fy  ; Value { Integral { [ CompY[CoefGeo * Tmax[{d a}] * N[] ]] ;
           In Domain ; Jacobian Vol ; Integration I1 ; } } }	   
     { Name F_vw ; Value {
         Integral { Type Global ; [ CoefGeo * 0.5 * nu[] * VirtualWork[{d a}] ];
           In ElementsOf[AirLayer, OnOneSideOf SkinDomainC_Moving];
           Jacobian Vol ; Integration I1 ; } } }


     { Name F_Lorentz ;
       Value { Integral { [
             CoefGeo * CrossProduct[-sigma[]*Dt[{a}], {d a}] ] ;
           In DomainC ;  Jacobian Vol; Integration I1; }
       }
     }

     //Only in Frequency domain  (steady-state)
     { Name F_cplx  ; Value { Integral { [  CoefGeo * Tmax_cplx[{d a}] * N[] ] ;
           In Domain ; Jacobian Vol ; Integration I1 ; } } }
     { Name Fmag_FD  ; Value { Term {  Type Global; [ Fmag_FD[] ] ; In DomainKin ; } } }

     { Name F_vw_cplx ; Value {
        Integral { Type Global ; [ CoefGeo * 0.5 * 0.5 * nu[] * VirtualWork[{d a}] ];
          In ElementsOf[AirLayer, OnOneSideOf SkinDomainC_Moving];
           Jacobian Vol ; Integration I1 ; } } }

     { Name Fmag_vw_FD  ; Value { Term {  Type Global; [ Fmag_vw_FD[] ] ; In DomainKin ; } } }

     { Name F_Lorentz_cplx ;
       Value { Integral { [
             CoefGeo * 0.5 * CrossProduct[Re[-sigma[]*Dt[{a}]], Re[{d a}]] +
             CoefGeo * 0.5 * CrossProduct[Im[-sigma[]*Dt[{a}]], Im[{d a}]]] ;
           In DomainC ;  Jacobian Vol; Integration I1; }
       }
     }
     { Name Florentz_FD  ; Value { Term {  Type Global; [ Florentz_FD[] ] ; In DomainKin ; } } }
     /*
     { Name F_Lorentz_cplx2 ; // Correct (same as previous)
       Value { Integral { [ CoefGeo * 0.5*Re[-sigma[]*Dt[{a}] /\ Conj[{d a}] ] ] ;
           In DomainC ;  Jacobian Vol; Integration I1; } } }
     */

     If(Flag_Cir)
       { Name Uc ; Value {
           Term { [ {Ub} ]  ; In DomainB ; }
           Term { [ {Uz} ]  ; In DomainZt_Cir ; }}}
       { Name Ic ; Value {
           Term { [ {Ib} ]  ; In DomainB ; }
           Term { [ {Iz} ]  ; In DomainZt_Cir ; }}}
     EndIf
   }
 }

 { Name Mechanical ; NameOfFormulation Mechanical ;
   PostQuantity {
     { Name U  ; Value { Term { [ {U} ]  ; In DomainKin ; } } } //Position
     { Name DeltaU  ; Value { Term { [ {U}-displacementY ]  ; In DomainKin ; } } } //Position change
     { Name V  ; Value { Term { [ {V} ]  ; In DomainKin ; } } } //Velocity
     { Name Fmag     ; Value { Term {  Type Global; [ Fmag[] ] ; In DomainKin ; } } }
     { Name Fmag_vw  ; Value { Term {  Type Global; [ Fmag_vw[] ] ; In DomainKin ; } } }
     { Name Ftot     ; Value { Term {  Type Global; [ Fmag[] + Fmg[] ] ; In DomainKin ; } } }
     { Name Florentz ; Value { Term {  Type Global; [ Florentz[] ] ; In DomainKin ; } } }
   }
 }

}


PostOperation MapMag_min UsingPost MagDyn_a_2D {
  Print[ az_nox, OnElementsOf Domain, File StrCat[ResDir,StrCat["tmp", ExtGmsh]], Format Gmsh,
    OverrideTimeStepValue 0, LastTimeStepOnly, SendToServer "No"] ;

  Print[ F[AirLayer], OnGlobal, Format Table, LastTimeStepOnly,
    File StrCat[ResDir, StrCat["Ftmp", ExtGnuplot]], Store 55] ;
///	 Print[ Fy[AirLayer], OnGlobal, Format Table, LastTimeStepOnly,
 //   File StrCat[ResDir, StrCat["Ftmp", ExtGnuplot]], Store 55] ;
  // Slow
   Print[ F_vw, OnRegion NodesOf[SkinDomainC_Moving], Format RegionValue, LastTimeStepOnly,
    File StrCat[ResDir, StrCat["Ftmp_vw", ExtGnuplot]], Store 56 ] ;
  Print[ F_Lorentz[DomainC_Moving], OnGlobal, Format Table, LastTimeStepOnly,
    File StrCat[ResDir, StrCat["Ftmp_lo", ExtGnuplot]], Store 57 ] ;
}

PostOperation MapMag_freq UsingPost MagDyn_a_2D {
  Print[ b, OnElementsOf Domain,  File StrCat[ResDir,StrCat["b", ExtGmsh]], Format Gmsh ] ;
  Print[ jz, OnElementsOf DomainC,  File StrCat[ResDir,StrCat["j", ExtGmsh]], Format Gmsh ] ;
  Print[ az, OnElementsOf Domain, File StrCat[ResDir,StrCat["az", ExtGmsh]], Format Gmsh ] ;
  Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
      "View[PostProcessing.NbViews-1].NbIso = 25;",
      "View[PostProcessing.NbViews-1].IntervalsType = 1;" // iso values
    ], File StrCat[ResDir,"az.opt"]];

  Print[ F_cplx[AirLayer], OnGlobal, Format Table, LastTimeStepOnly,
    File StrCat[ResDir, "Ftmp", ExtGnuplot], Store 65] ;
   Print[ Fmag_FD,    OnRegion DomainKin, Format Table, File > StrCat[ResDir, "Fmag_FD", ExtGnuplot],
   LastTimeStepOnly,SendToServer StrCat[output_mag,"41F_y MST [N]"], Color "LightYellow"] ;

  Print[ F_vw_cplx, OnRegion NodesOf[SkinDomainC_Moving], Format RegionValue, LastTimeStepOnly,
    File StrCat[ResDir, StrCat["Ftmp_vw", ExtGnuplot]], Store 67 ] ;
 Print[ Fmag_vw_FD, OnRegion DomainKin, Format Table, File > StrCat[ResDir,"Fmag_vw_cplx", ExtGnuplot],
    LastTimeStepOnly,SendToServer StrCat[output_mag,"42F_y VW [N]"], Color "LightYellow"] ;

  Print[ F_Lorentz_cplx[DomainC_Moving], OnGlobal, Format Table,
    File StrCat[ResDir, "F_lo_av.dat"], Store 66] ;
  Print[ Florentz_FD, OnRegion DomainKin, Format Table, File > StrCat[ResDir, "Florentz_FD", ExtGnuplot],
    LastTimeStepOnly, SendToServer StrCat[output_mag,"43F_y Lorentz law [N]"], Color "LightYellow"] ;
}


PostOperation MapMag_slow UsingPost MagDyn_a_2D {
  Print[ boundary,  OnElementsOf Dummy, File StrCat[ResDir, StrCat["bnd",ExtGmsh]], Format Gmsh,
    OverrideTimeStepValue step, LastTimeStepOnly ] ;
  Print[ b, OnElementsOf Domain,  File StrCat[ResDir,StrCat["b", ExtGmsh]], Format Gmsh,
    OverrideTimeStepValue step, LastTimeStepOnly ] ;
  Print[ jz, OnElementsOf DomainC,  File StrCat[ResDir,StrCat["j", ExtGmsh]], Format Gmsh,
    OverrideTimeStepValue step, LastTimeStepOnly ] ;
  Print[ az, OnElementsOf Domain, File StrCat[ResDir,StrCat["az", ExtGmsh]], Format Gmsh,
    OverrideTimeStepValue step, LastTimeStepOnly] ;
  Echo[Str["View[PostProcessing.NbViews-1].RangeType = 3;" ,// per timestep
      "View[PostProcessing.NbViews-1].NbIso = 25;",
      "View[PostProcessing.NbViews-1].IntervalsType = 1;" // iso values
    ], File StrCat[ResDir,"az.opt"]];
}

PostOperation MapMec UsingPost Mechanical {
  Print[ U, OnRegion DomainKin, File > StrCat[ResDir,StrCat["disp", ExtGnuplot]], Format Table,
    LastTimeStepOnly, SendToServer StrCat[output_mec, str_disp] ] ;
  Print[ DeltaU, OnRegion DomainKin, File StrCat[ResDir,StrCat["deltadisp", ExtGnuplot]], Format Table,
    LastTimeStepOnly, SendToServer "DeltaU"] ;
  Print[ V, OnRegion DomainKin, File > StrCat[ResDir,StrCat["velocity", ExtGnuplot]], Format Table,
    LastTimeStepOnly, SendToServer StrCat[output_mec, str_vel] ] ;

  Print[ Fmag,    OnRegion DomainKin, Format Table, File > StrCat[ResDir, StrCat["Fmag", ExtGnuplot]],
    LastTimeStepOnly,SendToServer StrCat[output_mag,"41F_y MST [N]"], Color "LightYellow"] ;
   Print[ Fmag_vw, OnRegion DomainKin, Format Table, File > StrCat[ResDir,StrCat["Fmag_vw", ExtGnuplot]],
    LastTimeStepOnly,SendToServer StrCat[output_mag,"42F_y VW [N]"], Color "LightYellow"] ;
  Print[ Florentz, OnRegion DomainKin, Format Table, File > StrCat[ResDir, StrCat["Florentz", ExtGnuplot]],
    LastTimeStepOnly,SendToServer StrCat[output_mag,"43F_y Lorentz law [N]"], Color "LightYellow"] ;
}



DefineConstant[
  // always solve this resolution (it's the only one provided)
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 0}
  // set some command-line options for getdp
  C_ = {"-solve -v 3 -v2 -bin", Name "GetDP/9ComputeCommand", Visible 0}
  // don't do the post-processing pass (some already included in Resolution)
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
];

If(Flag_AnalysisType==1)
  DefineConstant[
    // Some plots.... Far from handy !!! :-(
    // 1,2,3 -> Interval type
    xx_ = "2000", // Active X (second pair corresponds to X')
    yy_ = "0200", // Active Y (second pair corresponds to Y')
    null_="0000", // Not active for that graph

    time_ = { 0, Name "GetDP/Time", ReadOnly 1,
      Graph StrCat[ xx_,xx_, null_, xx_], Visible 1}

    // Levitation height as a function of time
    position_ = { displacementY, Name StrCat[output_mec, str_disp], ReadOnly 1,
      Graph StrCat[yy_,null_,xx_,null_], Visible 1}

    // Velocity as a function of time
    // Velocity as a function of height
    speed_    = { velocityY, Name StrCat[output_mec, str_vel], ReadOnly 1,
      Graph StrCat[null_, yy_, yy_, null_], Visible 1}

    // Force as a function of time
    force_ = { forceY, Name StrCat[output_mag,"41F_y MST [N]"], ReadOnly 1,
      Graph StrCat[null_,null_,null_,yy_], Visible 1, Highlight "LightYellow"}

    // WATCH OUT: this is kinda experimental... you could get strange results as the step for
    // computing the map is not well defined... NEEDS reworking
    // Some additional post-pro (slow, so to perform punctually)
    pos_maps = "c=DefineString['-pos MapMag_slow -bin -v 3',
    Name 'GetDP_NoAutoRun/9ComputeCommand', ReadOnly 1, NeverChanged 1, Visible 0];"
    one_run = "OnelabRun('GetDP_NoAutoRun',
    StrCat(Solver.Executable0, StrCat(' ', StrPrefix(General.FileName)))); Draw;"

    m_slow_postpro = { StrCat[pos_maps, one_run],
      Name StrCat[M_calc, "Maps: b, j, az"], Highlight "Wheat",
      Macro "GmshParseString", AutoCheck 0}
  ];
EndIf
If(Flag_AnalysisType==0)
  UndefineConstant[
    "GetDP/Time",
    StrCat[output_mec, str_vel],
    StrCat[output_mec, str_disp],
    StrCat[output_mag,"41F_y MST [N]"],
    StrCat[output_mag,"42F_y VW [N]"],
    StrCat[output_mag,"43F_y Lorentz law [N]"],
    StrCat[M_calc, "Maps: b, j, az"]
  ];
EndIf
