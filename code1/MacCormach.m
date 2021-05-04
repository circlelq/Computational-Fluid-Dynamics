function [] = MacCormach()
	global dt
	global dx
    global N
    global U
	global F
	global U2
	global U21

    for i = 2:N-1
        U21(i,:) = U(i,:) - dt/dx * (F(i,:)-F(i-1,:));
    end
    F = W2F(U2W(U21));
	for i = 2:N-1
        U2(i,:) = 0.5 * (U(i,:)+U21(i,:)) - 0.5 * dt/dx * (F(i+1,:)-F(i,:));
    end
    U = U2;
end

