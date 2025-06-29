function ID = findClosest(vector, value)
         D = abs(vector - value);
         Dmin = min(D);
         ID = find(D==Dmin);
end