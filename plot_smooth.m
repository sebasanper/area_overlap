z = load('final_speed_precise.dat');
clear x y fi
for d=0:0
for i=1:8
for j=1:10
k=8*(j-1)+i+d*80;
x(i,j)=z(k,2);
y(i,j)=z(k,3);
fi(i,j)=z(k,7);
end
end
surf(x,y,fi)
view(0,90)
shading interp
colormap jet
pause(0.1)
end