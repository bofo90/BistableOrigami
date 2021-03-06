function stop = outfun(x,optimValues,state)

stop = false;
 
switch state
   case 'init'
%        hold on
   case 'iter'
%        poop.x = ones(36,1)*optimValues.iteration;
       setGlobalx(x)
       % Concatenate current point and objective function
       % value with history. x must be a row vector.
%        history.fval = [history.fval; optimValues.fval];
%        history.x = [history.x; x];
       % Concatenate current search direction with 
       % searchdir.
%        searchdir = [searchdir;...
%                     optimValues.searchdirection'];
%        plot(x(1),x(2),'o');
       % Label points with iteration number.
       % Add .15 to x(1) to separate label from plotted 'o'
%        text(x(1)+.15,x(2),num2str(optimValues.iteration));
   case 'done'
%        hold off
   otherwise
end


function setGlobalx(val)
global poop
poop.x = [poop.x val];
%Turn the flag on only after making an iteration on the fmincon
poop.flag = 1;