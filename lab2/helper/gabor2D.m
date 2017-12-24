function h = gabor2D( U0,V0,a,b,Fs )
%Create a 2D gabor filter

N = floor(ceil((3*Fs/a)+1)/2);    %calculate size
[x,y] = meshgrid(-N:1:N,-N:1:N);  %make the grid

G = exp( - ( (a*x/Fs).^2 + (b*y/Fs).^2 ) ); %G(x,y)
C = cos( 2*pi*U0*x/Fs + 2*pi*V0*y/Fs );     %C(x,y)
S = sin( 2*pi*U0*x/Fs + 2*pi*V0*y/Fs );     %S(x,y)

h = G.*(C + 1i*S); %not normalized filter in space
h = h/norm(h,2);   %normalization of filter


end