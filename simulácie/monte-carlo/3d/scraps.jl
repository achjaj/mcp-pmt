function atan2(x, y)
    x > 0 && y >= 0 && return atan(y/x)
    x < 0 && return atan(y/x) + π
    x > 0 && y < 0 && return atan(y/x) + 2π
    x == 0 && y > 0 && return π/2
    x == 0 && y < 0 && return 3π/2 
end

cart2spher(vec) = [norm(vec), atan2(vec[1:2]...), acos(vec[3]/norm(vec))]