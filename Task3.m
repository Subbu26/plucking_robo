%%SUBRAMANIAN SIVANESAN 139495
syms L1 L2 L3 L4 L5 L6 Q1 Q2 Q3 Q4 Q5 Q6 m1 m2 m3 m4 m5 m6 Ixx1 Ixy1 Ixz1 Iyx1 Iyy1 Iyz1 Izx1 Izy1 Izz1 Ixx2 Ixy2 Ixz2 Iyx2 Iyy2 Iyz2 Izx2 Izy2 Izz2 Ixx3 Ixy3 Ixz3 Iyx3 Iyy3 Iyz3 Izx3 Izy3 Izz3 Ixx4 Ixy4 Ixz4 Iyx4 Iyy4 Iyz4 Izx4 Izy4 Izz4 Ixx5 Ixy5 Ixz5 Iyx5 Iyy5 Iyz5 Izx5 Izy5 Izz5 Ixx6 Ixy6 Ixz6 Iyx6 Iyy6 Iyz6 Izx6 Izy6 Izz6;
alpha=[0,-90,0,90,-90,90]; %%alpha value of links
a=[0,0,L2,L3,0,0];  %%lengQ of links
d=[0,L1,L4,L5,0,0];   %%offset
Q=[Q1,Q2,Q3,Q4,Q5,Q6];  %%Qeta represented by Q
for i=1:6
    switch i
        case 1
            T01= [cos(Q(1,i)),-sin(Q(1,i))*cosd(alpha(1,i)),sind(alpha(1,i))*sin(Q(1,i)),a(1,i)*cos(Q(1,i));sin(Q(1,i)),cos(Q(1,i)).*cosd(alpha(1,i)),-sind(alpha(1,i))*cos(Q(1,i)),sin(Q(1,i))*a(1,i);0,sind(alpha(1,i)),cosd(alpha(1,i)),d(1,i);0,0,0,1];
            
        case 2
            T12= [cos(Q(1,i)),-sin(Q(1,i))*cosd(alpha(1,i)),sind(alpha(1,i))*sin(Q(1,i)),a(1,i)*cos(Q(1,i));sin(Q(1,i)),cos(Q(1,i)).*cosd(alpha(1,i)),-sind(alpha(1,i))*cos(Q(1,i)),sin(Q(1,i))*a(1,i);0,sind(alpha(1,i)),cosd(alpha(1,i)),d(1,i);0,0,0,1];
            
        case 3
            T23= [cos(Q(1,i)),-sin(Q(1,i))*cosd(alpha(1,i)),sind(alpha(1,i))*sin(Q(1,i)),a(1,i)*cos(Q(1,i));sin(Q(1,i)),cos(Q(1,i)).*cosd(alpha(1,i)),-sind(alpha(1,i))*cos(Q(1,i)),sin(Q(1,i))*a(1,i);0,sind(alpha(1,i)),cosd(alpha(1,i)),d(1,i);0,0,0,1];
            
        case 4
            T34= [cos(Q(1,i)),-sin(Q(1,i))*cosd(alpha(1,i)),sind(alpha(1,i))*sin(Q(1,i)),a(1,i)*cos(Q(1,i));sin(Q(1,i)),cos(Q(1,i)).*cosd(alpha(1,i)),-sind(alpha(1,i))*cos(Q(1,i)),sin(Q(1,i))*a(1,i);0,sind(alpha(1,i)),cosd(alpha(1,i)),d(1,i);0,0,0,1];
            
         case 5
            T45= [cos(Q(1,i)),-sin(Q(1,i))*cosd(alpha(1,i)),sind(alpha(1,i))*sin(Q(1,i)),a(1,i)*cos(Q(1,i));sin(Q(1,i)),cos(Q(1,i)).*cosd(alpha(1,i)),-sind(alpha(1,i))*cos(Q(1,i)),sin(Q(1,i))*a(1,i);0,sind(alpha(1,i)),cosd(alpha(1,i)),d(1,i);0,0,0,1];
            
        case 6
            T56= [cos(Q(1,i)),-sin(Q(1,i))*cosd(alpha(1,i)),sind(alpha(1,i))*sin(Q(1,i)),a(1,i)*cos(Q(1,i));sin(Q(1,i)),cos(Q(1,i)).*cosd(alpha(1,i)),-sind(alpha(1,i))*cos(Q(1,i)),sin(Q(1,i))*a(1,i);0,sind(alpha(1,i)),cosd(alpha(1,i)),d(1,i);0,0,0,1]; 
    end
end

T01
T02 = T01*T12
T03 = T02*T23
T04 = T03*T34
T05 = T04*T45
T06 = T05*T56

T_matrix =[T01 T02 T03 T04 T05 T06 ];
Q_matrix =[Q1 Q2 Q3 Q4 Q5 Q6];  %%Qeta matrix

r1= transpose([-L1/2 0 0 1]); 
r2= transpose([-L2/2 0 0 1]); 
r3= transpose([-L3/2 0 0 1]); 
r4= transpose([-L4/2 0 0 1]); 
r5= transpose([-L5/2 0 0 1]);
r6= transpose([-L6/2 0 0 1]);

%%Uij matrices

   

   U11=diff(T01,Q1); U21=diff(T02,Q1); U31=diff(T03,Q1);
   U12=zeros(4);         U22=diff(T02,Q2); U32=diff(T03,Q2);
   U13=zeros(4);         U23=zeros(4);         U33=diff(T03,Q3);
   U14=zeros(4);         U24=zeros(4);         U34=zeros(4);         
   U15=zeros(4);         U25=zeros(4);         U35=zeros(4);         
   U16=zeros(4);         U26=zeros(4);         U36=zeros(4);
   
    U41=diff(T04,Q1); U51=diff(T05,Q1); U61=diff(T06,Q1);
   U42=diff(T04,Q2); U52=diff(T05,Q2); U62=diff(T06,Q2);
   U43=diff(T04,Q3); U53=diff(T05,Q3); U63=diff(T06,Q3);
   U44=diff(T04,Q4); U54=diff(T05,Q4); U64=diff(T06,Q4);
   U45=zeros(4);   U55=diff(T05,Q5); U65=diff(T06,Q5);
   U46=zeros(4);  U56=zeros(4);         U66=diff(T06,Q6);
   
   U111=diff(U11,Q1); U121=zeros(4)     ; U131=zeros(4);
   U112=zeros(4)     ; U122=zeros(4)     ; U132=zeros(4);
   U113=zeros(4)     ; U123=zeros(4)     ; U133=zeros(4);
   U114=zeros(4)     ; U124=zeros(4)     ; U134=zeros(4);
   U115=zeros(4)     ; U125=zeros(4)     ; U135=zeros(4);
   U116=zeros(4)     ; U126=zeros(4)     ; U136=zeros(4);
 
   U141=zeros(4)     ; U151=zeros(4)     ; U161=zeros(4);
   U142=zeros(4)     ; U152=zeros(4)     ; U162=zeros(4);
   U143=zeros(4)     ; U153=zeros(4)     ; U163=zeros(4);
   U144=zeros(4)     ; U154=zeros(4)     ; U164=zeros(4);
   U145=zeros(4)     ; U155=zeros(4)     ; U165=zeros(4);
   U146=zeros(4)     ; U156=zeros(4)     ; U166=zeros(4);
 
%%Uijk matrices

%%U2jk matrices
   U211=diff(U21,Q1); U221=diff(U22,Q1); U231=zeros(4);
   U212=diff(U21,Q2); U222=diff(U22,Q2); U232=zeros(4);
   U213=zeros(4)     ; U223=zeros(4)     ; U233=zeros(4);
   U214=zeros(4)     ; U224=zeros(4)     ; U234=zeros(4);
   U215=zeros(4)     ; U225=zeros(4)     ; U235=zeros(4);
   U216=zeros(4)     ; U226=zeros(4)     ; U236=zeros(4);
  
   
   U241=zeros(4)     ; U251=zeros(4)     ; U261=zeros(4);
   U242=zeros(4)     ; U252=zeros(4)     ; U262=zeros(4);
   U243=zeros(4)     ; U253=zeros(4)     ; U263=zeros(4);
   U244=zeros(4)     ; U254=zeros(4)     ; U264=zeros(4);
   U245=zeros(4)     ; U255=zeros(4)     ; U265=zeros(4);
   U246=zeros(4)     ; U256=zeros(4)     ; U266=zeros(4);
%U3jk Matrixes
   U311=diff(U31,Q1); U321=diff(U32,Q1); U331=diff(U33,Q1);
   U312=diff(U31,Q2); U322=diff(U32,Q2); U332=diff(U33,Q2);
   U313=diff(U31,Q3); U323=diff(U32,Q3); U333=diff(U33,Q3);
   U314=zeros(4)     ; U324=zeros(4)     ; U334=zeros(4);
   U315=zeros(4)     ; U325=zeros(4)     ; U335=zeros(4);
   U316=zeros(4)     ; U326=zeros(4)     ; U336=zeros(4);
   U341=zeros(4)     ; U351=zeros(4)     ; U361=zeros(4);
   U342=zeros(4)     ; U352=zeros(4)     ; U362=zeros(4);
   U343=zeros(4)     ; U353=zeros(4)     ; U363=zeros(4);
   U344=zeros(4)     ; U354=zeros(4)     ; U364=zeros(4);
   U345=zeros(4)     ; U355=zeros(4)     ; U365=zeros(4);
   U346=zeros(4)     ; U356=zeros(4)     ; U366=zeros(4);
   
   U411=diff(U41,Q1); U421=diff(U42,Q1); U431=diff(U43,Q1);
   U412=diff(U41,Q2); U422=diff(U42,Q2); U432=diff(U43,Q2);
   U413=diff(U41,Q3); U423=diff(U42,Q3); U433=diff(U43,Q3);
   U414=diff(U41,Q4); U424=diff(U42,Q4); U434=diff(U43,Q4);
   U415=zeros(4);      U425=zeros(4);      U435=zeros(4);
   U416=zeros(4);      U426=zeros(4);      U436=zeros(4);
   
   U441=diff(U44,Q1); U451=zeros(4);      U461=zeros(4);      
   U442=diff(U44,Q2); U452=zeros(4);      U462=zeros(4);      
   U443=diff(U44,Q3); U453=zeros(4);      U463=zeros(4);      
   U444=diff(U44,Q4); U454=zeros(4);      U464=zeros(4);      
   U445=zeros(4);      U455=zeros(4);      U465=zeros(4);      
   U446=zeros(4);      U456=zeros(4);      U466=zeros(4); 
   
   U511=diff(U51,Q1); U521=diff(U52,Q1); U531=diff(U53,Q1); 
   U512=diff(U51,Q2); U522=diff(U52,Q2); U532=diff(U53,Q2); 
   U513=diff(U51,Q3); U523=diff(U52,Q3); U533=diff(U53,Q3); 
   U514=diff(U51,Q4); U524=diff(U52,Q4); U534=diff(U53,Q4); 
   U515=diff(U51,Q5); U525=diff(U52,Q5); U535=diff(U53,Q5); 
   U516=zeros(4);      U526=zeros(4);      U536=zeros(4);           
   
   U541=diff(U54,Q1); U551=diff(U55,Q1); U561=zeros(4);  
   U542=diff(U54,Q2); U552=diff(U55,Q2); U562=zeros(4);
   U543=diff(U54,Q3); U553=diff(U55,Q3); U563=zeros(4);      
   U544=diff(U54,Q4); U554=diff(U55,Q4); U564=zeros(4);      
   U545=diff(U54,Q5); U555=diff(U55,Q5); U565=zeros(4);      
   U546=zeros(4);      U556=zeros(4);      U566=zeros(4);
   
   U611=diff(U61,Q1); U621=diff(U62,Q1); U631=diff(U63,Q1); 
   U612=diff(U61,Q2); U622=diff(U62,Q2); U632=diff(U63,Q2); 
   U613=diff(U61,Q3); U623=diff(U62,Q3); U633=diff(U63,Q3); 
   U614=diff(U61,Q4); U624=diff(U62,Q4); U634=diff(U63,Q4); 
   U615=diff(U61,Q5); U625=diff(U62,Q5); U635=diff(U63,Q5); 
   U616=diff(U61,Q6); U626=diff(U62,Q6); U636=diff(U63,Q6);
   U641=diff(U64,Q1); U651=diff(U65,Q1); U661=diff(U66,Q1); 
   U642=diff(U64,Q2); U652=diff(U65,Q2); U662=diff(U66,Q2);
   U643=diff(U64,Q3); U653=diff(U65,Q3); U663=diff(U66,Q3);
   U644=diff(U64,Q4); U654=diff(U65,Q4); U664=diff(U66,Q4);
   U645=diff(U64,Q5); U655=diff(U65,Q5); U665=diff(U66,Q5);
   U646=diff(U64,Q6); U656=diff(U65,Q6); U666=diff(U66,Q6);
 %I matrixes
 
    I1=[Ixx1 Ixy1 Ixz1;
        Iyx1 Iyy1 Iyz1;
        Izx1 Izy1 Izz1];
   
    I2=[Ixx2 Ixy2 Ixz2;
        Iyx2 Iyy2 Iyz2;
        Izx2 Izy2 Izz2];
    
    I3=[Ixx3 Ixy3 Ixz3;
        Iyx3 Iyy3 Iyz3;
        Izx3 Izy3 Izz3];
    
    I4=[Ixx4 Ixy4 Ixz4;
        Iyx4 Iyy4 Iyz4;
        Izx4 Izy4 Izz4];
    
    I5=[Ixx5 Ixy5 Ixz5;
        Iyx5 Iyy5 Iyz5;
        Izx5 Izy5 Izz5];
    
    I6=[Ixx6 Ixy6 Ixz6;
        Iyx6 Iyy6 Iyz6;
        Izx6 Izy6 Izz6]
    
  %J matrix
  
    I=zeros(4);
  for i=1:6
      I=eval(['I' num2str(i)]);
      m=eval(['m' num2str(i)]);
      r=eval(['r' num2str(i)]);
      eval(['j11' '=((-I(1,1)+I(2,2)+I(3,3))/2)']);
      eval(['j12' '=I(1,2)']);
      eval(['j13' '=I(1,3)']);
      eval(['j14' '=m*r(1)']);
      
      eval(['j21' '=I(1,2)']);
      eval(['j22' '=((I(1,1)-I(2,2)+I(3,3))/2)']);
      eval(['j23' '=I(2,3)']);
      eval(['j24' '=m*r(2)']);
      
      eval(['j31' '=I(1,3)']);
      eval(['j32' '=I(2,3)']);
      eval(['j33' '=((I(1,1)+I(2,2)-I(3,3))/2)']);
      eval(['j34' '=m*r(3)']);
      
      eval(['j41' '=m*r(1)']);
      eval(['j42' '=m*r(2)']);
      eval(['j43' '=m*r(3)']);
      eval(['j44' '=m']);
      
      J=[j11 j12 j13 j14;j21 j22 j23 j24; j31 j32 j33 j34;j41 j42 j43 j44];
      eval(['J' num2str(i) '=J']);
  end
  
   %D matrix
 Uaux=zeros(4);
 for i=1:6
   for j=1:6
       m=max([i j]);
       x=1;
       for k=m:6
       Uaux=eval(['U' num2str(k) num2str(i)]);
       Ud=Uaux';
       A=eval(['U' num2str(k) num2str(j)])*eval(['J' num2str(k)])*Ud;
       x=x+trace(A);
       end
       eval(['d' num2str(i) num2str(j) '=x']);
   end
 end

D=[d11 d12 d13 d14 d15 d16 ;
   d21 d22 d23 d24 d25 d26 ;
   d31 d32 d33 d34 d35 d36 ;
   d41 d33 d43 d44 d45 d46 ;
   d51 d52 d53 d54 d55 d56 ;
   d61 d62 d63 d64 d65 d66 ];

%Hikm matrix
  for i=1:6
   for k=1:6
     for m=1:6
       j=max([i m k]);
       x=0;
       for l=j:6
       Uaux=eval(['U' num2str(j) num2str(i)]);
       Uh=Uaux';
       x=x+trace(eval(['U' num2str(j) num2str(k) num2str(m)])*eval(['J' num2str(j)])*Uh);
       end
       eval(['h' num2str(i) num2str(k) num2str(m) '=x']);
     end 
   end
  end
  
    %Coriolis column matrix
 syms dQ1 dQ2 dQ3 dQ4 dQ5 dQ6
 
  for i=1:6
      y=0;
      for k=1:6 
        y=y+x;
        x=0;  
          for m=1:6
          y=x+eval(['h' num2str(i) num2str(k) num2str(m)])*eval(['dQ' num2str(k)])*eval(['dQ' num2str(m)]);
          end    
      end
      eval(['h' num2str(i) '=y']);
  end
  
   H=[h1;h2;h3;h4;h5;h6];
   
     %Gravity column matrix
  
  g=[0 0 -9.81 0];
  
    for i=1:6
      x=0;
      for j=i:6
      x=x+(-eval(['m' num2str(j)])*g*eval(['U' num2str(j) num2str(i)])*eval(['r' num2str(j)]));  
      end
      eval(['c' num2str(i) '=x']);
    end
  
      C=[c1;c2;c3;c4;c5;c6];
  
  %Final dynamic model
  
  %angular acceleration matrix
  
  syms ddQ1 ddQ2 ddQ3 ddQ4 ddQ5 ddQ6
  ddQ=[ddQ1;ddQ2;ddQ3;ddQ4;ddQ5;ddQ6]; 
  %T=[t1;t2;t3;t4;t5;t6]; %Matrix wiQ torQues of each joint
  
  T=D*ddQ+H+C
  