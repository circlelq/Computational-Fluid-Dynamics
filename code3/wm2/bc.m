%
%  bc.m
%  code3
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/6/11.
%

function [W] = bc(W)
% set boundary conditions
%   
    global Nx Ny
	d = (Nx*0.2);
	h = (Ny*0.2);

	% Entry condition
	for j = 1:Ny
		W(1, j, :) = [1.4, 3, 0, 1];
	end

	% Export condition
	for j = 1:Ny
		W(Nx, j, :) = W(Nx-2, j, :);
	end

	% Top condition
	for i = 2:Nx-1
		W(i, Ny, :) = W(i, Ny-2, :);
		W(i, Ny, 3) = -W(i, Ny-2, 3);
	end

	% Bottom condition
	for i = 2:d
		W(i, 1, :) = W(i, 3, :);
		W(i, 1, 3) = -W(i, 3, 3);
	end

	for i = d+1:Nx-1
		W(i, h, :) = W(i, h+2, :);
		W(i, h, 3) = -W(i, h+2, 3);
	end


	for j = 1:h
		W(d+1, j, :) = W(d-1, j, :);
		W(d+1, j, 2) = -W(d-1, j, 2);
	end

	% corner condition
	W(d+1, h, 1) = (W(d+1, h+2, 1)+W(d-1, h, 1))/2;
	W(d+1, h, 2) = (W(d+1, h+2, 2)-W(d-1, h, 2))/2;
	W(d+1, h, 3) = (-W(d+1, h+2, 3)+W(d-1, h, 3))/2;
	W(d+1, h, 4) = (W(d+1, h+2, 4)+W(d-1, h, 4))/2;






end

