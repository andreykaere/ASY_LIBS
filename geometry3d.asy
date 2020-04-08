import three;
import graph3;
import solids;

import geometry;

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

/*
path3 circle3(triple A, triple B, triple C) {
    return 
}
*/
path project3(path3 p) {
    return path(p, project3);
}


bool collinear3(triple A, triple B) {
    return dot(A,B) == 0;
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





