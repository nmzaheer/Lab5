Include "t28_data.geo";

cen = newp; Point(newp)= {0,0,0,lr0};
//Inner Coil points
ptIn[]+= newp; Point(newp) = {ic_radius-width_ic/2,0,0,lci};
ptIn[]+= newp; Point(newp) = {ic_radius+width_ic/2,0,0,lci};
ptIn[]+= newp; Point(newp) = {ic_radius+width_ic/2,-height_coil,0,lci};
ptIn[]+= newp; Point(newp) = {ic_radius-width_ic/2,-height_coil,0,lci};

lcIn[]+= newl; Line(newl) = {ptIn[0],ptIn[1]};
lcIn[]+= newl; Line(newl) = {ptIn[1],ptIn[2]};
lcIn[]+= newl; Line(newl) = {ptIn[2],ptIn[3]};
lcIn[]+= newl; Line(newl) = {ptIn[3],ptIn[0]};

Line Loop(newll) = {lcIn[]};
surfInCoil[]+= news; Plane Surface(news) = {newll-1};

//Outer Coil points
ptOut[]+= newp; Point(newp) = {oc_radius-width_oc/2,0,0,lci};
ptOut[]+= newp; Point(newp) = {oc_radius+width_oc/2,0,0,lci};
ptOut[]+= newp; Point(newp) = {oc_radius+width_oc/2,-height_coil,0,lci};
ptOut[]+= newp; Point(newp) = {oc_radius-width_oc/2,-height_coil,0,lci};

lcOut[]+= newl; Line(newl) = {ptOut[0],ptOut[1]}; 
lcOut[]+= newl; Line(newl) = {ptOut[1],ptOut[2]}; 
lcOut[]+= newl; Line(newl) = {ptOut[2],ptOut[3]}; 
lcOut[]+= newl; Line(newl) = {ptOut[3],ptOut[0]};

Line Loop(newll) = {lcOut[]};
surfOutCoil[]+= news; Plane Surface(news) = {newll-1};
//Plate points
ptPlate[]+= newp; Point(newp) = {0,initial_height,0,lcp};
ptPlate[]+= newp; Point(newp) = {0,initial_height+plate_thick,0,lcp};
ptPlate[]+= newp; Point(newp) = {plate_radius,initial_height+plate_thick,0,lcp};
ptPlate[]+= newp; Point(newp) = {plate_radius,initial_height,0,lcp};

lp[]+= newl; Line(newl) = {ptPlate[0],ptPlate[1]};
lp[]+= newl; Line(newl) = {ptPlate[1],ptPlate[2]};
lp[]+= newl; Line(newl) = {ptPlate[2],ptPlate[3]};
lp[]+= newl; Line(newl) = {ptPlate[3],ptPlate[0]};

Line Loop(newll) = {lp[]};
surfPlate[]+= news; Plane Surface(news) = {newll-1};

// Air shell
ptShell[]+= newp; Point(newp) = {0,Rro,0,lr0};
ptShell[]+= newp; Point(newp) = {0,Rri,0,lr0};
ptShell[]+= newp; Point(newp) = {0,-Rri,0,lr0};
ptShell[]+= newp; Point(newp) = {0,-Rro,0,lr0};
ptShell[]+= newp; Point(newp) = {Rri,0,0,lr0};
ptShell[]+= newp; Point(newp) = {Rro,0,0,lr0};

lsi[]+= newl; Circle(newl) = {ptShell[1],cen,ptShell[4]};
lsi[]+= newl; Circle(newl) = {ptShell[4],cen,ptShell[2]};

lso[]+= newl; Circle(newl) = {ptShell[0],cen,ptShell[5]};
lso[]+= newl; Circle(newl) = {ptShell[5],cen,ptShell[3]};

//AirLayer for force computation
ptAl[]+= newp; Point(newp) = {0,initial_height-width_al,0,lal};
ptAl[]+= newp; Point(newp) = {plate_radius+width_al,initial_height-width_al,0,lal};
ptAl[]+= newp; Point(newp) = {plate_radius+width_al,initial_height+plate_thick+width_al,0,lal};
ptAl[]+= newp; Point(newp) = {0,initial_height+plate_thick+width_al,0,lal};

lnal[]+= newl; Line(newl) = {ptAl[0],ptAl[1]};
lnal[]+= newl; Line(newl) = {ptAl[1],ptAl[2]};
lnal[]+= newl; Line(newl) = {ptAl[2],ptAl[3]};
lnal[]+= newl; Line(newl) = {ptAl[3], ptPlate[1]};
lnal[]+= newl; Line(newl) = {ptPlate[0],ptAl[0]};
lnal[]+= newl; Line(newl) = {ptAl[1],ptPlate[3]};
lnal[]+= newl; Line(newl) = {ptAl[2],ptPlate[2]};

//Air Layer Srface Up
Line Loop(newll) = {lp[1],-lnal[6],lnal[{2,3}]};
surfAirLayUp[]+= news; Plane Surface(news) = {newll-1};

//Air Layer Surface Right
Line Loop(newll) = {lp[2],-lnal[5],lnal[{1,6}]};
surfAirLayRight[]+= news; Plane Surface(news) = {newll-1};

//Air Layer Surface Bottom
Line Loop(newll) = {lnal[4],lnal[0],lnal[5],lp[3]};
surfAirLayBot[]+= news; Plane Surface(news) = {newll-1};

//Axisymmetry line
ls[]+=newl; Line(newl) = {ptShell[0],ptShell[1]};
ls[]+=newl; Line(newl) = {ptShell[1],ptAl[3]};
ls[]+=newl; Line(newl) = {ptAl[0],cen};
ls[]+=newl; Line(newl) = {cen,ptShell[2]};
ls[]+=newl; Line(newl) = {ptShell[2],ptShell[3]};

Line Loop(newll) = {ls[0],lsi[],ls[4],-lso[]};
surfAirShell[]+=news; Plane Surface(news) = {newll-1}; //Outer Air shell

Line Loop(newll) = {ls[1],-lnal[{2,1,0}],ls[{2,3}],-lsi[],lcIn[],lcOut[]};
surfInShell[]+=news; Plane Surface(news) = {newll-1}; // Inner Air shell

Physical Surface(INNER_COIL) = surfInCoil[];
Physical Surface(OUTER_COIL) = surfOutCoil[];
Physical Surface(PLATE) = surfPlate[];
Physical Surface(AIR_SHELL) = surfAirShell[];
Physical Surface(INNER_SHELL) = surfInShell[];
Physical Surface(AIR_LAY_UP) = surfAirLayUp[];
Physical Surface(AIR_LAY_RIGHT) = surfAirLayRight[];
Physical Surface(AIR_LAY_DOWN) = surfAirLayBot[];
Physical Line(OUTER_BND) = lso[];
Physical Line(AXIS_SYM) = {ls[{0,1}],lnal[3],-lp[0],lnal[4],ls[{2,3,4}]};
Physical Line(SKIN_DOMAIN) = {lp[{1,2,3}]};
