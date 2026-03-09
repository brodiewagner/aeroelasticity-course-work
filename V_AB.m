function [VAB] = V_AB(A,B,C)
   
    r0 = [B(1) - A(1), B(2) - A(2), B(3) - A(3)] ;
    r1 = [C(1) - A(1), C(2) - A(2), C(3) - A(3)] ;
    r2 = [C(1) - B(1), C(2) - B(2), C(3) - B(3)] ;

    % Use dot product to find the squared magnitude of the cross product
    n = cross(r1, r2) / dot(cross(r1,r2), cross(r1,r2)) ;

    % Use norm() to find scalar magnitude
    VAB = (1/(4*pi)) .* n * (dot(r0, (r1/norm(r1)) - (r2/norm(r2)))) ;
end