import three;
import graph3;
import solids;

import geometry;

filltype dotfilltype = Fill;


real length3(triple A, triple B=O) {
    return sqrt((A.x - B.x)^2 + (A.y - B.y)^2 + (A.z - B.z)^2);
}


triple[] get_basis(projection P = currentprojection) {
	triple Zp = unit(P.camera);
    triple Xp = unit(cross(Z, Zp));
    triple Yp = unit(cross(Zp, Xp));
    triple[] basis = {Xp, Yp, Zp};
    
    return basis;
}



//pair project3(triple A, projection P = currentprojection) {
pair project3(triple A) {
    projection P = currentprojection;
    triple[] basis = get_basis(P);
    triple Xp = basis[0];
    triple Yp = basis[1];
    triple Zp = basis[2];

    real[][] transform_matrix = {{Xp.x, Yp.x, Zp.x},
                                 {Xp.y, Yp.y, Zp.y},
                                 {Xp.z, Yp.z, Zp.z}
                                };
    real[] A_matrix = {A.x, A.y, A.z};
    
    real[] Ap = solve(transform_matrix, A_matrix);
    
    return (Ap[0], Ap[1]);    
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
               filltype filltype_ = NoFill,
               margin margin = NoMargin, marker marker = nomarker) {
    path3 w = Circle(B,radius, normal(A--B--C));
    triple P = intersectionpoint(w,B--A); //line3 BA
    triple Q = intersectionpoint(w,B--C); //line3 BC
    path3 arc = Arc(B,P,Q);
    path3 g = B--P--arc--B--cycle;
    draw(pic, L, arc, p, arrow=arrow3);
    fill(project3(g), filltype_.fillpen);
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
triple project2(pair A, path3 p) {
    projection P = currentprojection;
    triple[] basis = get_basis(P);
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

bool collinear3(triple A, triple B) {
    return dot(A,B) == 0;
}


path3 Circle(triple C, triple A, triple normal=Z){
    return Circle(C,length3(C,A),normal);
}



struct line3 {

    triple A,B;
    
    path3 line;
    path3 line_extended;
    
    triple[] inits;   

    void operator init(triple A, triple B) {
        
    }

    bool on_line (triple P) {
        return collinear3(P-A, P-B);
    }

    line3 xline3 (real x) {
        line3 l;    

        return l;
    }
    
}

struct plane3 {
    real A,B,C,D;
    surface surface;
    surface surface_extended;
    triple [] inits;
    path3 contour;        

    void operator init(triple A, triple B, triple C) {
        
        //this.inits = {A,B,C};
    }

    bool in_plane (triple P) {
        return this.A*P.x + this.B*P.y + this.C*P.z + this.D == 0;
    }


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
    
    triple Pa3 = invert(Pa, normal(A--B--C)); 
    triple Pb3 = invert(Pb, normal(A--B--C)); 
    triple Pc3 = invert(Pc, normal(A--B--C)); 

    draw(Pa3--Pb3--Pc3);
}



triple midpoint3(triple A, triple B) {
    return (A+B)/2;
}


triple foot(triple T, triple A, triple B, triple C) {
    return midpoint3(T, reflect(A,B,C)*T);
}


/*
triple selectpoint3(triple A, triple B, ) {

}
*/


triple foot(triple A, plane3 a) {
    triple[] inits = a.inits;
    return midpoint3(A, reflect(inits[0], inits[1], inits[2])*A);
}

/*
triple foot(triple A, triple P, triple Q) {
    real d = abs(cross(unit(P-Q)));
    line3.line_extended 
}
*/

triple foot(triple A, line3 l) {
    triple P = l.inits[0];
    triple Q = l.inits[1];
    
    real d = abs(cross(unit(P-Q), unit(A-P))*length3(A, P));
    
    return intersectionpoint(l.line_extended, 
                             Circle(A, d, normal(A--P--Q))); 
}



transform3 orthogonalproject(plane3 p) {
	triple[] inits = p.inits;	
	path3 g = inits[0]--inits[1]--inits[2]--cycle;
	return planeproject(g);
}


path3 incircle3(triple A, triple B, triple C) {return nullpath3;}

path3 circle3(triple A, triple B, triple C){return nullpath3;}

/*
void markangle3{}



path3 line3  {} ...

...(surface) or path3  plane(triple A, triple B, triple C, sizze ....) {}





foot of  point ... on plane/line
*/





