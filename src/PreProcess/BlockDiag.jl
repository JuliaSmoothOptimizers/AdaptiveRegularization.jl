n=5
vD1 = rand(n)
vD2 = zeros(n-1)
vD2[1] = 1.0
vD2[3] = 1.0

D = spdiagm((vD1,vD2,vD2),[0 -1 1]) # don't forget inquire returns D^-1
#DiagD, Q = eig(full(D))  

vQ1 = ones(vD1)
vQ2 = zeros(vD2)
veig = copy(vD1)

i=1;
while i<n
    if vD2[i] == 0.0
        #veig[i] = vD1[i]
        #vQ1[i] = 1.0
        #vQ2[i] = 0.0
        i += 1
    else
        mA = [vD1[i] vD2[i];vD2[i] vD1[i+1]]
        DiagmA, Qma = eig(mA)
        veig[i] = DiagmA[1]
        vQ1[i] = Qma[1,1]
        vQ2[i] = Qma[1,2]
        vQ1[i+1] = Qma[2,2]
        veig[i+1] = DiagmA[2]
        i += 2
    end  
end

Q = spdiagm((vQ1,vQ2,vQ2),[0,-1,1])
