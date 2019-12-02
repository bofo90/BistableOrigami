function angles = anglesDesign4Vertex(opt)

switch opt.vertexType
    case {'non'}
        angles = opt.angDesign;
    case{'random'}
        ang = randAng();
        angles = [0, ang(1), sum(ang(1:2)), sum(ang)]*pi/180;
        
    case{'X1'}
        ang = randAngOrdered(1);
        ang = ang([1 2 4 3]);        
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'X2'}
        ang = randAngOrdered(2);
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'Y1'}
        ang = randAngOrdered(1);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'Y2'}
        ang = randAngOrdered(2);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;    
	case{'Z1'}
        ang = randAngOrdered(1);
        ang = ang([1 3 2 4]);        
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'Z2'}
        ang = randAngOrdered(2);
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;     
        
    case{'2NS'}
        ang = randAng2Equal('s');
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'2OS'}
        ang = randAng2Equal('s');
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
        
    case{'2NL'}
        ang = randAng2Equal('l');
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
    case{'2OL'}
        ang = randAng2Equal('l');
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
        
    case{'2OM1'}
        ang = randAng2Equal('m', 1);
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
    case{'2OM2'}
        ang = randAng2Equal('m', 2);
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
    case{'2NM1'}
        ang = randAng2Equal('m', 1);
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
    case{'2NM2'}
        ang = randAng2Equal('m', 2);
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;
        
    case{'GFF'}
        ang = randAngSumPi();
        ang = ang([1 2 4 3]);        
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'CY'}
        ang = randAngSumPi();
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;
    case{'CZ'}
        ang = randAngSumPi();
        ang = ang([1 3 2 4]);        
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;
        
    case{'CFF'}
        ang = randAngSumPiEqual();
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
    case{'2C'}
        ang = randAngSumPiEqual();
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
        
    case{'2OFF'}
        ang = randAngSumPiOver2();
        ang = ang([1 2 4 3]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
    case{'2NC'}
        ang = randAngSumPiOver2();
        ang = ang([1 3 2 4]);
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180; 
        
    case{'3S'}
        ang = randAng3Equal('s');
        ang = ang([1 2 4 3]);        
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
    case{'3L'}
        ang = randAng3Equal('l');
        ang = ang([1 2 4 3]);        
        angles = [0, ang(1), sum(ang(1:2)), sum(ang(1:3))]*pi/180;  
        
    case{'2CFF'}
        angles = [0, 90, 180, 270]*pi/180; 
        
    otherwise
        error('\n----------\nThe type of 4-Vertex is non-existent or specified incorrectly\n----------\n',[])

end   

function ang = randAng()
while 1
    ang = rand(1,3)*150+15;
    ang(4) = 360 - sum(ang);
    if ~(sum(ang>15) < 4 || sum(ang<165) < 4)
        break;
    end
end

function ang = randAngOrdered(subtype)
while 1
    ang = rand(1,3)*150+15;
    ang(4) = 360 - sum(ang);
    ang = sort(ang);
    if ~(sum(ang>15) < 4 || sum(ang<165) < 4)
        if (subtype == 1) && (ang(1)+ang(4) < ang(2)+ang(3))
            break;
        end
        if (subtype == 2) && (ang(1)+ang(4) > ang(2)+ang(3))
            break;
        end
    end
end

function ang = randAng2Equal(equal, subtype)
while 1
    ang = rand(1,2)*150+15;
    ang(3) = ang(1);
    ang(4) = 360 - sum(ang);
    ang = sort(ang);
    if ~(sum(ang>15) < 4 || sum(ang<165) < 4)
        if equal == 's' && ang(1) == ang(2)
            break;
        end
        if equal == 'l' && ang(3) == ang(4)
            break;
        end
        if equal == 'm' && ang(2) == ang(3)
            if (subtype == 1) && (ang(1)+ang(4) < ang(2)+ang(3))
                break;
            end
            if (subtype == 2) && (ang(1)+ang(4) > ang(2)+ang(3))
                break;
            end
        end
    end
end

function ang = randAngSumPi()
ang = rand(1,2)*150+15;
ang(3) = 180 - ang(2);
ang(4) = 180 - ang(1);
ang = sort(ang);

function ang = randAngSumPiEqual()
ang = rand(1,1)*150+15;
ang(2) = ang(1);
ang(3) = 180 - ang(2);
ang(4) = ang(3);
ang = sort(ang);

function ang = randAngSumPiOver2()
ang = rand(1,1)*150+15;
ang(2) = 90;
ang(3) = 90;
ang(4) = 180 - ang(1);
ang = sort(ang);

function ang = randAng3Equal(equal)
while 1
    ang = rand(1,1)*150+15;
    ang(2) = ang(1);
    ang(3) = ang(1);
    ang(4) = 360 - sum(ang);
    ang = sort(ang);
    if ~(sum(ang>15) < 4 || sum(ang<165) < 4)
        if equal == 's' && ang(1) == ang(2)
            break;
        end
        if equal == 'l' && ang(3) == ang(4)
            break;
        end
    end
end