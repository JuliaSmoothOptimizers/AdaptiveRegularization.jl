f = :foo
offset = 2

@eval begin
    function $f(x)
        println("This is one arg ", $f)
        x + $offset
    end
    
    function $f(x,y)
        println("This is two args ", $f)
        return $f(x+y)
    end
    
end

@eval function  $f(x,y,z)
println("This is three args ", $f)
end
