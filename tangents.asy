import geometry;

line[] out_tangents(circle w1, circle w2) {
    var r1 = w1.r;
    var r2 = w2.r;
    point o1 = w1.C;
    point o2 = w2.C;
        
    var dy = o1.y - o2.y;
    var dx = o1.x - o2.x;   
    var l  = length(segment(o1, o2));
    
    if (r1 == r2) {
        circle I = circle(o2, sqrt(l*l + r1*r1));
        point A = intersectionpoints(w1, I)[0];
        point B = intersectionpoints(w1, I)[1];
        
        line[] ans = {tangent(w1, A), tangent(w1, B)};
        return ans;
    }
    
    var  dM = l + (l * r1/r2) / (1 - r1/r2);
    point M = (o2.x + dx * dM/l, o2.y + dy * dM/l) ;
    return tangents(w1, M);
}

line [] in_tangents(circle w1, circle w2) {
    var r1 = w1.r;
    var r2 = w2.r;
    point o1 = w1.C;
    point o2 = w2.C;
        
    var dy = o1.y - o2.y;
    var dx = o1.x - o2.x;   
    var l  = length(segment(o1, o2));

    var dS = (l * r1/r2) / (1 + r1/r2); //distance between o1 and 
                                        //intersection point of tangents

    point M = intersectionpoints(circle(o1, dS), segment(o1, o2))[0];
    //point M = (o1.x - dx * dS/l , o1.y - dy * dS/l);
    return tangents(w1, M);
}


