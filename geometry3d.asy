import three;
import graph3;
import solids;

import geometry;

filltype dotfilltype = Fill;



struct basis3 {
    triple x,y,z;

    void operator init(triple x,triple y,triple z) {
        this.x = unit(x);
        this.y = unit(y);
        this.z = unit(z);
    }
}


struct vector3 {
    real x,y,z;
    
    void operator init(triple V) {
        this.x = V.x;
        this.y = V.y;
        this.z = V.z;
    }   
    
    void operator init(triple A, triple B) {
        this.operator init(B-A);
    }   

}


triple operator ecast(vector3 v) {
    return (v.x,v.y,v.z);
}

vector3 operator *(vector3 v, real k) {
    return vector3((k*v.x, k*v.y,k*v.z));
}

triple operator +(triple A, vector3 v) {
    return A + (triple) v;
}

triple operator -(triple A, vector3 v) {
    return A - (triple) v;
}

vector3 cross(vector3 v, vector3 u) {
    return vector3(cross((triple) v, (triple) u));
}

vector3 dot(vector3 v, vector3 u) {
    return dot((triple) v, (triple) u);
}

real length3(vector3 v) {
    return (v.x)^2 + (v.y)^2 + (v.z)^2;
}

bool collinear3(vector3 v, vector3 u) {
    return length3(cross(v,u)) == 0;
}



struct curve3 {
    path[] front;
    path[] back;
}


struct object3 {
    string type;
    surface surface;
    curve3 curve;
    bool drawn;

    void operator init(string type,
                       surface surface, 
                       curve3 curve) {
        this.curve = curve;
        this.surface = surface;
        this.type = type;
    }
};


object3[] OBJECTS;


struct ray3 {
    vector3 vec;
    triple A;
    curve3 curve;
    path3 ray;  

    void operator init(triple A, vector3 v, real k=1000) {
        this.vec = v;
        this.A = A;
        this.ray = A--(A+v*k);
    }


    void operator init(triple A, triple B) {
        this.operator init(A, vector3(A,B));
    }

}


struct segment3 {}


struct line3 {

	triple A,B;
    triple extendA, extendB;
	path3 line;
	path3 line_extended;
    vector3 vec;
	curve3 curve;   

    void operator init(triple A, triple B, real k=1000) {
        this.vec = vector3(A,B); //vector AB
        this.A = A;
        this.B = B;
        this.extendB = B + vec * k;
        this.extendA = A - vec * k;
        this.line_extended = extendA--extendB;    
        this.line = A--B;
    }
    
    void operator init(vector3 v, triple A) {
        this.operator init(A, A+v);
    } 
   
    //replace with operator
    /*
    line3 xline3 (real x) {
        line3 l;    

        return l;
    }
    */

}




struct plane3 {
    real A,B,C,D;
    surface surface;
    surface surface_extended;
    triple[] inits;
    curve3 curve;
    real length;
    real width;    
    vector3 normal;    
    triple center;
    path3 contour;

    void operator init(triple C, vector3 normal, 
                        real length=50, real width=50, real angle=0) {
        this.normal = normal;        
        this.center = C;        
        //this.inits = {A,B,C};
        this.length = length;
        this.width = width;
        
        vector3 camera = vector3(currentprojection.camera);
    
        if (collinear3(camera, normal)) {}
        else {}
 
        vector3 oX = unit(cross(camera,normal)) * width/2;
        vector3 oY = unit(cross(oX, normal)) * length/2;
        
        path3 contour = (C+oX+oY)--(C+oX-oY)--(C-oX-oY)--(C-oX+oY)--cycle;
        contour = rotate(angle, normal) * contour;
     
        this.contour = contour;
        //rotate by the angle "angle" about v

        curve.front = project3(contour);
        
        
        
    }



    void operator init(triple A, triple B, triple C, 
                        real length=50, real width=50, real angle=0) { //angle in degrees, if angle = 0 it means that one side of plane
                                                                      // is parallel to axis oX
        this.normal = vector3(cross(A-B,A-C));        
        this.center = intersectionpoint((A--midpoint3(B,C)),
                                        (B--midpoint3(A,C)));       
        
        this.inits.push(A);
        this.inits.push(B);
        this.inits.push(C);

        this.length = length;
        this.width = width;
        



        object3 object = object3("surface",surface,curve);            
        OBJECTS.push(object);          
            
    }
    
    


    void operator init(line3 a, triple A) {
        if (A @ a) abort("point lies on the line, can't create plane,
                          passing through them");
        this.operator init(A, a.A, a.B);       
    }


}


struct plane_line {}




//-----------------------------------------------------------------------
//functions



real distance3(triple A, triple B) {
    return sqrt((A.x - B.x)^2 + (A.y - B.y)^2 + (A.z - B.z)^2);
}

triple midpoint3(triple A, triple B) {
    return (A+B)/2;
}


basis3 getBasis(projection P = currentprojection) {
    triple Zp = P.camera;
    triple Xp = cross(P.up, Zp);
    triple Yp = cross(Zp, Xp);
    
    return basis3(Xp,Yp,Zp);
}


triple calcCoordsInBasis(basis3 basis, triple A) {
    triple Xp = basis.x;
    triple Yp = basis.y;
    triple Zp = basis.z;
    
    triple Ap = (Xp.x*A.x + Yp.x*A.y + Zp.x*A.z,
                 Xp.y*A.x + Yp.y*A.y + Zp.y*A.z,
                 Xp.z*A.x + Yp.z*A.y + Zp.z*A.z);
    return Ap;

}


triple changeBasis(basis3 basis1, basis3 basis2, triple A) {
    triple X2 = basis2.x;
    triple Y2 = basis2.y;
    triple Z2 = basis2.z;
    
    triple Xp = calcCoordsInBasis(basis1, X2);
    triple Yp = calcCoordsInBasis(basis1, Y2);
    triple Zp = calcCoordsInBasis(basis1, Z2);

    real[][] transform_matrix =  {{Xp.x, Yp.x, Zp.x},
                                  {Xp.y, Yp.y, Zp.y},
                                  {Xp.z, Yp.z, Zp.z}
                                 };
    
    real[] A_matrix = {A.x, A.y, A.z};
    
    real[] Ap = solve(transform_matrix, A_matrix);

    return (Ap[0],Ap[1],Ap[2]);
    
}


//pair project3(triple A, projection P = currentprojection) {
pair project3(triple A) {
    //projection P = currentprojection;
    //triple Ap = changeBasis(basis3(X,Y,Z),getBasis(P),A);    
    //return (Ap.x, Ap.y);
    return project(A);
} 

path project3(path3 p) {
    return path(p, project3);
}



//copied from https://tex.stackexchange.com/questions/321674/how-to-make-empty-point-in-asymptote which is copied from three_surface.asy
void opendot(picture pic=currentpicture, triple v, material p=currentpen,
         light light=nolight, string name="", render render=defaultrender)
{
  pen q=(pen) p;
  pen fillpen = light.background;
  if (invisible(fillpen)) fillpen = currentlight.background;
  if (invisible(fillpen)) fillpen = white;
  real size=0.5*linewidth(dotsize(q)+q);
  pic.add(new void(frame f, transform3 t, picture pic, projection P) {
      triple V=t*v;
      assert(!is3D(), "opendot() not supported unless settings.prc == false and settings.render != 0");
      if(pic != null)
        dot(pic,project(V,P.t),filltype=FillDraw(fillpen=fillpen, drawpen=q));
    },true);
  triple R=size*(1,1,1);
  pic.addBox(v,v,-R,R);
}

void opendot(picture pic=currentpicture, Label L, triple v, align align=NoAlign,
             string format=defaultformat, material p=currentpen,
             light light=nolight, string name="", render render=defaultrender)
{
  Label L=L.copy();
  if(L.s == "") {
    if(format == "") format=defaultformat;
    L.s="("+format(format,v.x)+","+format(format,v.y)+","+
      format(format,v.z)+")";
  }
  L.align(align,E);
  L.p((pen) p);
  opendot(pic,v,p,light,name,render);
  label(pic,L,v,render);
}
//

void markangle3(picture pic = currentpicture,
               Label L = "", int n = 1, real radius = 0, real space = 0,
               explicit triple A, explicit triple B, explicit triple C,
                 pair align = dir(1),
               arrowbar3 arrow3 = None, pen p = currentpen,
               filltype filltype = NoFill,
               margin margin = NoMargin, marker marker = nomarker) {
    path3 w = Circle(B,radius, normal(A--B--C));
    triple P = intersectionpoint(w,B--A); //line3 BA
    triple Q = intersectionpoint(w,B--C); //line3 BC
    path3 arc = Arc(B,P,Q);
    path3 g = B--P--arc--B--cycle;
    draw(pic, L, arc, p, arrow=arrow3);
    fill(project3(g), filltype.fillpen);
}
/*
void markangle3(picture pic = currentpicture,
               Label L = "", int n = 1, real radius = 0, real space = 0,
               explicit line3 l1, explicit line3 l2, explicit pair align = dir(1),
               arrowbar arrow = None, pen p = currentpen,
               filltype filltype = NoFill,
               margin margin = NoMargin, marker marker = nomarker) {}

*/
/*
triple invert3(pair A, path3 p) {
    projection P = currentprojection;
    triple[] basis = getBasis(P);
    triple Xp = basis[0];
    triple Yp = basis[1];
    triple Zp = basis[2];

    real[][] transform_matrix = {{Xp.x, Yp.x, Zp.x},
                                 {Xp.y, Yp.y, Zp.y},
                                 {Xp.z, Yp.z, Zp.z}
                                };
    real[] A_matrix = {A.x, A.y, 0};
    
    real[] Ap = solve(transform_matrix, A_matrix);
    
    return (Ap[0], Ap[1]);

}
*/

/*
don not need to uncomment

void dot3(picture pic=currentpicture, Label L, triple A, align align=NoAlign,
         string format=defaultformat, pen p=currentpen, filltype filltype) {

    dot(pic,L,project3(A),align,format,p,filltype);
}

void draw3(picture pic=currentpicture, Label L="", path3 g, 
          align align=NoAlign, pen p=currentpen) {

    draw(pic,L,project3(g),align,p);

}

void draw3(picture pic=currentpicture, Label L="", guide3 g, 
          align align=NoAlign, pen p=currentpen) {

    draw(pic,L,project3(g),align,p);

}

void draw3(picture pic=currentpicture, Label L="", path3[] G, 
          align align=NoAlign, pen p=currentpen) {

    for(path3 g:G) {
        draw3(pic,L,g,align,p);
    }
}
*/



/*
path3 circle3(triple A, triple B, triple C) {
    return 
}
*/





bool collinear3(triple A, triple B, triple C) {
    return collinear3(vector3(A,B), vector3(A,C));
}


path3 Circle(triple C, triple A, triple normal=Z){
    return Circle(C,distance3(C,A),normal);
}





line3 parallel(line3 a, triple A) {
    return line3(a.vec, A);
}



triple[] intersectionpoints(line3 a, surface s) {
    return intersectionpoints(a.line_extended,s);
}


bool is_intersecting(line3 a, surface s) {
    triple[] u = intersectionpoints(a, s);
    return u.length > 0;
}

/*
triple intersectionpoint(line3 a, plane3 s) {
    //return intersectionpoints(a.line_extended,s);
    if (a.A @ s) {return a.A}
    
    //matrix A = solve();...
}
*/

bool is_intersecting(line3 a, plane3 s) {
    triple u = intersectionpoint(a, s);
    return u != nullpath3 ;
}



bool operator @(triple P, line3 a) {
	return collinear3(a.A,a.B,P);
}


line3 invertpoint(pair A, projection P=currentprojection) {
    vector3 vec = vector3(P.camera);
    
    triple[] basis = getBasis(P);
    triple Xp = basis[0];
    triple Yp = basis[1];
    triple Zp = basis[2];
    triple A = (A.x,A.y,1);
    /*
    real[][] transform_matrix = {{Xp.x, Yp.x, Zp.x},
                                 {Xp.y, Yp.y, Zp.y},
                                 {Xp.z, Yp.z, Zp.z}
                                };
    */
    //real[] A_matrix = {A.x, A.y, 1};
    
    triple Ap = (Xp.x*A.x + Yp.x*A.y + Zp.x*A.z,
                 Xp.y*A.x + Yp.y*A.y + Zp.y*A.z,
                 Xp.z*A.x + Yp.z*A.y + Zp.z*A.z);
    
    //real[] Ap = solve(transform_matrix, A_matrix);
    //triple Ap = (Ap[0],Ap[1],Ap[2]);
    //triple Ap = invert(A, (0,0,1)); //just get one (random) point on the line
    
    return line3(Ap,Ap+vec);
} //get a line from pair



line3 invertpoint(triple A, projection P=currentprojection) {
    return line3(A, A+P.camera);
}








//stolen from module geometry.asy
void Drawline(picture pic = currentpicture,Label L = "", line3 a,  pair P, bool dirP = true, pair Q, bool dirQ = true,
                      align align = NoAlign, pen p = currentpen,
                      arrowbar arrow = None,
                      Label legend = "", marker marker = nomarker,
                      pathModifier pathModifier = NoModifier)
{/* Add the two parameters 'dirP' and 'dirQ' to the native routine
    'drawline' of the module 'math'.
    Segment [PQ] will be prolonged in direction of P if 'dirP = true', in
    direction of Q if 'dirQ = true'.
    If 'dirP = dirQ = true', the behavior is that of the native 'drawline'.
    Add all the other parameters of 'Draw'.*/
    
  pic.add(new void (frame f, transform t, transform T, pair m, pair M) {
      picture opic;
      // Reduce the bounds by the size of the pen.
      m -= min(p) - (linemargin(), linemargin()); M -= max(p) + (linemargin(), linemargin());

      // Calculate the points and direction vector in the transformed space.
      t = t * T;
      pair z = t * P;
      pair q = t * Q;
      pair v = q - z;
      // path g;
      pair ptp, ptq;
      real cp = dirP ? 1:0;
      real cq = dirQ ? 1:0;
      // Handle horizontal and vertical lines.
      if(v.x == 0) {
        if(m.x <= z.x && z.x <= M.x)
          if (dot(v, m - z) < 0) {
            ptp = (z.x, z.y + cp * (m.y - z.y));
            ptq = (z.x, q.y + cq * (M.y - q.y));
          } else {
            ptq = (z.x, q.y + cq * (m.y - q.y));
            ptp = (z.x, z.y + cp * (M.y - z.y));
          }
      } else if(v.y == 0) {
        if (dot(v, m - z) < 0) {
          ptp = (z.x + cp * (m.x - z.x), z.y);
          ptq = (q.x + cq * (M.x - q.x), z.y);
        } else {
          ptq = (q.x + cq * (m.x - q.x), z.y);
          ptp = (z.x + cp * (M.x - z.x), z.y);
        }
      } else {
        // Calculate the maximum and minimum t values allowed for the
        // parametric equation z + t * v
        real mx = (m.x - z.x)/v.x, Mx = (M.x - z.x)/v.x;
        real my = (m.y - z.y)/v.y, My = (M.y - z.y)/v.y;
        real tmin = max(v.x > 0 ? mx : Mx, v.y > 0 ? my : My);
        real tmax = min(v.x > 0 ? Mx : mx, v.y > 0 ? My : my);
        pair pmin = z + tmin * v;
        pair pmax = z + tmax * v;
        if(tmin <= tmax) {
          ptp = z + cp * tmin * v;
          ptq = z + (cq == 0 ? v:tmax * v);
        }
      }
      path g = ptp--ptq;
      if (length(g)>0)
        {
          if(L.s != "") {
            Label lL = L.copy();
            if(L.defaultposition) lL.position(Relative(.9));
            lL.p(p);
            lL.out(opic, g);
          }
          g = pathModifier(g);
          if(linetype(p).length == 0){
            pair m = midpoint(g);
            pen tp;
            tp = dirP ? p : addpenline(p);
            draw(opic, pathModifier(m--ptp), tp);
            tp = dirQ ? p : addpenline(p);
            draw(opic, pathModifier(m--ptq), tp);
          } else {
            draw(opic, g, p);
          }

    
          marker.markroutine(opic, marker.f, g);
          arrow(opic, g, p, NoMargin);
          add(f, opic.fit());
        }
    });
}


void draw(picture pic=currentpicture, line3 a, Label L = "", bool dirA=true, 
            bool dirB=true, bool inf=true, pen p=currentpen) {
    

    if (inf) {
        Drawline(pic, L, a, project3(a.A), dirA, project3(a.B), dirB, p);
    }	
    else {
        draw(pic,L,a.line,p);
        //write();    
    }
}
  



plane3 perpendicular(triple A, line3 a, 
					 real length=50, real width=50) {
    return plane3(A,a.vec,length,width); 
}




//A*x + B*y + C*z + D = 0
triple getpointXY(real x, real y, plane3 a) {
    real z;
    if (a.C == 0) {z = 0;}
    else {z = (a.A*x + a.B*y + a.D)/a.C;}
    
    return (x,y,z);
}

triple getpointXZ(real x, real z, plane3 a) {
    real y;
    if (a.B == 0) {y = 0;}
    else {y = (a.A*x + a.C*z + a.D)/a.B;}
    
    return (x,y,z);
}

triple getpointYZ(real y, real z, plane3 a) {
    real x;
    if (a.A == 0) {x = 0;}
    else {x = (a.B*y + a.C*z + a.D)/a.A;}
    
    return (x,y,z);
}


bool operator @(triple Q, plane3 a) { //, bool inf=true) {
    if (inf) {
        return a.A*Q.x + a.B*Q.y + a.C*Q.z + a.D == 0;
    }
	
    //return a.A*Q.x + a.B*Q.y + a.C*Q.z + a.D == 0 && 
    //       is_intersecting(invertpoint(Q), a);
}



bool is_skew(line3 a, line3 b) {
    return !(b.A @ plane3(a, b.B));
}

bool is_parallel(line3 a, line3 b) {
    return collinear3(a.vec, b.vec);
}

bool is_intersecting(line3 a, line3 b) {
    return !(is_skew(a,b)) && !(is_parallel(a,b));
}

triple intersectionpoint(line3 a, line3 b) {
    if (!is_intersecting(a,b)) {
        abort("lines do not intersect");
    }
    
    return intersectionpoint(a.line_extended,b.line_extended);
}

triple invert3 (pair A, line3 a) {
    return intersectionpoint(invertpoint(A),a);
}

triple intersectionpoint(line3 a, plane3 s, bool inf=true) {
    if (inf) {
        return intersectionpoints(a, s.surface_extended)[0];
    }
    
    return intersectionpoints(a,s.surface)[0];
}


triple getpointX(real x, line3 a) {
    if (a.vec.x == 0) { 
        abort("Line is parallel to the plane ZY, 
                so relust point can't be unique"); // handle special case
    }
    plane3 planeZYp = plane3((x,0,0),
                             (x,1,0),
                             (x,0,1)); //plane, which is parallel to planeZY and intersecting axis oX at x
    
    return intersectionpoint(a, planeZYp);
        
    
}

triple getpointX(real y, line3 a) {
    if (a.vec.y == 0) { 
        abort("Line is parallel to the plane XZ, 
                so relust point can't be unique"); // handle special case
    }
    
    plane3 planeXZp = plane3((0,y,1),
                             (1,y,0),
                             (0,y,1)); //plane, which is parallel to planeXZ and intersecting axis oY at y
    
    return intersectionpoint(a, planeXZp);
        
    
}


triple getpointZ(real z, line3 a) {
    if (a.vec.z == 0) { 
        abort("Line is parallel to the plane XY, 
                so relust point can't be unique"); // handle special case
    }
    
    plane3 planeXYp = plane3((1,0,z),
                             (0,1,z),
                             (0,0,z)); //plane, which is parallel to planeXY and intersecting axis oZ at z
    
    return intersectionpoint(a, planeXYp);
        
    
}









//void markrightangle3 {}
/*
void markrightangle3(triple ptI, triple ptO, triple ptJ, real size=0.5, pen p=currentpen) {
   triple imI=ptO+size*unit(ptI-ptO),
          imJ=ptO+size*unit(ptJ-ptO),
          imK=imI+imJ-ptO;
   draw(imI--imK--imJ,p); 
}
*/

void markrightangle3(triple A, triple B, triple C, real n=5, pen p=currentpen) {
    real s = n/10;

    pair Ap = project3(A);     
    pair Bp = project3(B);     
    pair Cp = project3(C);
    
    //dot(Ap);
    //dot(Bp);
    //dot(Cp);

    path w = circle(Bp, s);
    //draw(w);
    pair Pa = intersectionpoint((Bp--Ap),w);
    pair Pc = intersectionpoint((Bp--Cp),w);
    
    //dot(Pa);
    //dot(Pc);
    
    pair Pb = reflect(line(Pa, Pc)) * Bp;
    
    //replace 
    triple Pa3 = invert(Pa, normal(A--B--C)); 
    triple Pb3 = invert(Pb, normal(A--B--C)); 
    triple Pc3 = invert(Pc, normal(A--B--C)); 
    //

    draw(Pa3--Pb3--Pc3);
}





/*
triple foot3(triple T, triple A, triple B, triple C) {
    return midpoint3(T, reflect(A,B,C)*T);
}
*/

/*
triple selectpoint3(triple A, triple B, ) {

}
*/


triple foot3(triple A, plane3 a) {
    triple[] inits = a.inits;
    return midpoint3(A, reflect(inits[0], inits[1], inits[2])*A);
}

/*
triple foot(triple A, triple P, triple Q) {
    real d = abs(cross(unit(P-Q)));
    line3.line_extended 
}
*/

triple foot3(triple A, line3 l) {
    triple P = l.A;
    triple Q = l.B;
    
    real d = abs(cross(unit(P-Q), unit(A-P))*distance3(A, P));
    
    return intersectionpoint(l.line_extended, 
                             Circle(A, d, normal(A--P--Q))); 
}



transform3 orthogonalproject(plane3 p) {
	triple[] inits = p.inits;	
	path3 g = inits[0]--inits[1]--inits[2]--cycle;
	return planeproject(g);
}


//circle3
circle3 incircle3(triple A, triple B, triple C) {return nullpath3;}


//init
circle3 circle3(triple A, triple B, triple C){return nullpath3;}

/*
void markangle3{}



path3 line3  {} ...

...(surface) or path3  plane(triple A, triple B, triple C, sizze ....) {}





foot of  point ... on plane/line
*/


void drawCurve(picture pic=currentpicture, curve3 curve, 
                pen frontpen=currentpen, pen backpen=currentpen+dashed) {
    for (path front : curve.front) {
        draw(pic, front, frontpen);
    }
    
    for (path back : curve.back) {
        draw(pic, back, backpen);
    }
}


void draw3(picture pic=currentpicture, line3 a, Label L = "", 
            bool dirB=true, bool inf=true, pen p=currentpen){

        object3 object = object3("surface",surface,curve);            
        OBJECTS.push(object);          
}

void withGeometry3d(void main()) {
    main();
    projection P = currentprojection;
    if (P.camera == P.up) {
        abort("camera's vector coincides with up-vector");
    } 

    add2dFrame();

    drawAllObjects();
}


/*
void create3D() {
	exitfcn currentexitfunction=atexit();
	atexit(exitfunction);
}
*/
