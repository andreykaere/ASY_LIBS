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
triple invert3(pair A, path3 p) {
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
    triple extendA, extendB;
	path3 line;
	path3 line_extended;
    triple vec;

	triple[] inits = {A,B};   

    void operator init(triple A, triple B, real k=1.3) {
        this.vec = B-A; //vector AB
        this.A = A;
        this.B = B;
        this.extendB = B+vec*k;
        this.extendA = A-vec*k;
        this.line_extended = extendA--extendB;    
        this.line = A--B;
    }
    
    //replace with operator
    bool on_line (triple P) {
        return collinear3(P-A, P-B);
    }
    /*
    line3 xline3 (real x) {
        line3 l;    

        return l;
    }
    */

}


line3 invertpoint(pair A, projection P=currentprojection) {
    triple vec = P.camera;
    
    triple[] basis = get_basis(P);
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


triple intersectionpoint(line3 a, line3 b) {
    return intersectionpoint(a.line_extended,b.line_extended);
}


triple invert3 (pair A, line3 a) {
    return intersectionpoint(invertpoint(A),a);
}


void Drawline(picture pic = currentpicture, Label L = "", pair P, bool dirP = true, pair Q, bool dirQ = true,
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


void draw(picture pic=currentpicture, line3 a, bool dirA=true, 
            bool dirB=true, bool inf=true, pen p=currentpen) {
    

    if (inf) {
        Drawline(project3(a.A),dirA,project3(a.B),dirB,p);
    }	
    else {
        draw(a.A--a.B,p);    
    }
}
  /*  
  pic.add(new void (frame f, transform t, transform T, pair m, pair M) {
      // Reduce the bounds by the size of the pen.
      m -= min(p); M -= max(p);

      triple P,Q;
      // Calculate the points and direction vector in the transformed space.
      t=t*T;
      pair z=t*project3(a.A);
      pair v=t*project3(a.B)-z;
      // Handle horizontal and vertical lines.
      if(v.x == 0) {
        if(m.x <= z.x && z.x <= M.x)
          P = invert3((z.x,m.y), a);
          Q = invert3((z.x,M.y), a);
          draw(f,P--Q,p);
        
      } else if(v.y == 0) {
        if(m.y <= z.y && z.y <= M.y)
          P = invert3((m.x,z.y), a);
          Q = invert3((M.x,z.y), a);
          draw(f,P--Q,p);
      } else {
        // Calculate the maximum and minimum t values allowed for the
        // parametric equation z + t*v
        real mx=(m.x-z.x)/v.x, Mx=(M.x-z.x)/v.x;
        real my=(m.y-z.y)/v.y, My=(M.y-z.y)/v.y;
        real tmin=max(v.x > 0 ? mx : Mx, v.y > 0 ? my : My);
        real tmax=min(v.x > 0 ? Mx : mx, v.y > 0 ? My : my);
        if(tmin <= tmax)
          dot(inverse(t)*(z+tmin*v));
          write(inverse(t)*(z+tmin*v));
          write("min");
          write(inverse(t)*(z+tmax*v));
          write("max");
          //dot((z+tmin*v));
          dot((3.86575335008085,-0.0117484232145937));
          dot((-0.176338164471332,5.38018313243324));
          dot((-0.0117484232145937,5.16062931569675));
          //dot(f,z+tmax*v);
		  //dot(invertpoint(inverse(t)*(z+tmin*v)).A,blue);
		  //write(invertpoint(z+tmin*v).A);
		  //dot(f,invertpoint(z+tmax*v).B,blue);
*/
/*
          P = invert3(inverse(t)*(z+tmin*v), a);
          Q = invert3(inverse(t)*(z+tmin*v), a);
          draw(f,P--Q,p);
  */
/*
    }

       
    },true);

}
*/

    //size(pic,0,0);
    //picture pic1;
    //picture pic2;
    //transform t = inverse(pic.calculateTransform());
  /*  
path bbox(picture pic=currentpicture,
            real xmargin=0, real ymargin=xmargin,
            pen p=currentpen, filltype filltype=NoFill)
 {
   frame f=pic.fit(max(pic.xsize-2*xmargin,0),max(pic.ysize-2*ymargin,0));
   return box(f,xmargin,ymargin,p,filltype,above=false);
}
*/
/*
	pic.add(new void (frame f, transform t, transform T, pair m, pair M) {
		picture opic;
		// Reduce the bounds by the size of the pen.
		m -= min(p); 
		M -= max(p);
		path border = box(m,M);

		path line_2d = project3(a.line_extended);
        draw((t*T)*border);
		//pair[] u = intersectionpoints(border,line_2d);

		//dot(opic, u[0]);
		//dot(opic, u[1]);
		add(f, opic.fit());
    });    
  */ 
/*
    if (inf) {
        pic.add(new void(frame f, transform3 t, picture pic1, projection P) {
			
			//picture pic1;
            frame f1 = bbox(pic);
            path b = box(f1);
            draw(f1,b);
            //path border = box(f1,bbox());//box((min(bbox(pic))-min(p)),(max(bbox(pic))-max(p)));
            //path border = box(m,M);
			//draw(border, blue);
			//draw(pic1,a.line_extended, green);
			draw(project3(a.line_extended), red);
			//dot(a.extendA);
			//draw(a.line,red);
			//write(a.A);
			//write(a.B);
			
			//pair[] u = intersectionpoints(project3(a.line_extended), bbox());

			//pair P = u[0];
			//pair Q = u[1];
			
			//dot(pic,P);
			
			//triple U = intersectionpoint(invertpoint(P),a);
			//triple V = intersectionpoint(invertpoint(Q),a);
			
			//dot(U);
			//draw(pic1, U--V, p);
            //add(pic, pic1.fit());
  */              
            //add(pic, f1);
            /*

            write(b); 
            draw(f1, (-106.897256650345,-106.851187223828)--(106.897256650345,-106.851187223828)--(106.897256650345,113.987345113472)--(-106.897256650345,113.987345113472)--cycle);

            add(pic, f1);
            write(is3D(f1));
            draw(box(f1),purple);
            write(min(inverse(pic.calculateTransform())*f1),max(f1));           
*/            

            // path border2 = inverse(pic1.calculateTransform())*box((min(bbox(pic))-min(p)),(max(bbox(pic))-max(p)));
            //draw(border2,orange);   
    // },true);
    //} else {
        //draw(pic1,a.A--a.B);
    //}
    //add(pic.fit3());
    //add(pic, pic1.fit());
    //add(pic, pic2.fit());
    
    //draw(project3(a.line_extended));
//}


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
    
    //replace 
    triple Pa3 = invert(Pa, normal(A--B--C)); 
    triple Pb3 = invert(Pb, normal(A--B--C)); 
    triple Pc3 = invert(Pc, normal(A--B--C)); 
    //

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
/*
triple foot(triple A, line3 l) {
    triple P = l.inits[0];
    triple Q = l.inits[1];
    
    real d = abs(cross(unit(P-Q), unit(A-P))*length3(A, P));
    
    return intersectionpoint(l.line_extended, 
                             Circle(A, d, normal(A--P--Q))); 
}
*/


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



void add_2d_frame() {
    picture pic1;
	draw(pic1,box(bbox()), invisible);

	add(pic1.fit());
}



void with_geometry3d(void main()) {
    main();
    add_2d_frame();
}


/*
void create3D() {
	exitfcn currentexitfunction=atexit();
	atexit(exitfunction);
}
*/
