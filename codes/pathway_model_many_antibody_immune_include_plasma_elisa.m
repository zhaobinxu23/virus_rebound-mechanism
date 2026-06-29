function F=pathway_model_many_antibody_immune_include_plasma_elisa(t,y,para,para_new,para_new_1)

 y = max(0,y);
 F = zeros(1402,1);


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j,1) = 0;
      
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+100,1) = -5*para_new(i)*y(10*(i-1)+j+100)*y(1402)+para_new_1(j)*y(10*(i-1)+j+900);
            
        
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+200,1) = 0;
        
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+300,1) =  -para_new(i)*y(10*(i-1)+j+300)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1100);
               

    end
end


%% c1-c4

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+400,1) =  0;

    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+500,1) = 0;
        F_E_2(i,j) = 0;
        
    end
end


for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+600,1) =  0;
        F_E_3(i,j) = 0;

    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+700,1) = 0;
        F_E_4(i,j) = 0;
    end
end
%% c1'-c4'

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+800,1) = 0;
        F_V_1(i,j) = 0;    
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+900,1) = 5*para_new(i)*y(10*(i-1)+j+100)*y(1402)-para_new_1(j)*y(10*(i-1)+j+900);
        F_V_2(i,j) = -5*para_new(i)*y(10*(i-1)+j+100)*y(1402)+para_new_1(j)*y(10*(i-1)+j+900);         
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1000,1) = 0;
        F_V_3(i,j) = 0;      
    end
end

for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1100,1) = para_new(i)*y(10*(i-1)+j+300)*y(1402)-para_new_1(j)*y(10*(i-1)+j+1100);
        F_V_4(i,j) = -para_new(i)*y(10*(i-1)+j+300)*y(1402)+para_new_1(j)*y(10*(i-1)+j+1100);      
    end
end

%% plasma_M
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1200,1) = 0;
             
    end
end

%% plasma_G
for i = 1:10
    for j = 1:10
        F(10*(i-1)+j+1300,1) = 0;
    end
end





F(1401,1) = 0;
F(1402,1) = sum(sum(F_V_1))+sum(sum(F_V_2))+sum(sum(F_V_3))+sum(sum(F_V_4));




end

