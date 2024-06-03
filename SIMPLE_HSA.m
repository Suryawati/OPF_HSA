clc
clear all;

%=========================parameter HSA ===================================
NVAR=2;         % number of variables
NG=6;           % number of ineguality constraints
NH=0;           % number of eguality constraints
MaxItr= 100 ;    % maximum number of iterations
HMS=40;         % harmony memory size
HMCR=0.9;       % harmony consideration rate  0< HMCR <1
PARmin=0.4;     % minumum pitch adjusting rate
PARmax=0.9;     % maximum pitch adjusting rate
bwmin=0.0001;   % minumum bandwidth
bwmax=1.0;      % maxiumum bandwidth
PVB=[0 5;
     0 5];
% /**** Initiate Matrix ****/
HM=zeros(HMS,NVAR);
NCHV=zeros(1,NVAR);
BestGen=zeros(1,NVAR);
fitness=zeros(1,HMS);
BW=zeros(1,NVAR);
gx=zeros(1,NG);

%==========================================================================
%======================main program========================================

        
        for i=1:HMS
            for j=1:NVAR
                %HM(i,j)=randval(PVB(j,1),PVB(j,2));
                HM(i,j)=rand(1)*(PVB(j,2)-PVB(j,1))+PVB(j,1);
            end
            %fitness(i) = Fitness(HM(i,:));
            fitness(i)=100*(HM(i,2)-HM(i,1))^2 + (HM(i,1)-1)^2;
        end
        
        currentIteration  = 0;
        
        while currentIteration < MaxItr
            
            PAR=(PARmax-PARmin)/(MaxItr)*currentIteration+PARmin;
            coef=log(bwmin/bwmax)/MaxItr;
            for pp =1:NVAR
                BW(pp)=bwmax*exp(coef*currentIteration);
            end
            
            % improvise a new harmony vector (NCHV) ==============================
        for i=1:HMS
            for j =1:NVAR
                ran = rand(1);                
                if( ran < HMCR ) % memory consideration
                    index = round(rand(1)*(HMS-1)+1);
                    NCHV(i,j) = HM(index,j);
                    pvbRan=rand(1);
                    if( pvbRan < PAR) % pitch adjusting
                        pvbRan1 = rand(1);
                        result = NCHV(i,j);
                        if( pvbRan1 < 0.5)
                            result =result+  rand(1) * BW(j);
                            if( result < PVB(j,2))
                                NCHV(i,j) = result;
                            end                        
                    
                        else                        
                            result =result- rand(1) * BW(j);                                                
                            if( result > PVB(j,1))
                                NCHV(i,j) = result;
                            end
                        end
                    end                        
                                    
                else
                    NCHV(i,j) = rand(1)*(PVB(j,2)-PVB(j,1))+PVB(j,1); % random selection
                end                          
            end
        end
        %calculate new fitness
            newFitness=zeros(HMS,1);
            for i=1:HMS
            newFitness(i) = 100*(NCHV(i,2)-NCHV(i,1))^2 + (NCHV(i,1)-1)^2;
            end
            
                      
            if(currentIteration==0)  
                [BestFit, BestIndex]=min(fitness);
                [WorstFit]=max(fitness);            
            end
            
        for i=1:HMS
        if (newFitness(i) < WorstFit)
            if( newFitness(i) < fitness(i) )
                HM(i,:)=NCHV(i,:);
                fitness(i)=newFitness(i);                
            else
               HM(i,:)=HM(BestIndex,:);
               fitness(i)=BestFit;
            end
        end           
             
        end
        
      [BestFit, BestIndex]=min(fitness);
      [WorstFit, WorstIndex]=max(fitness);  
      BestGen=HM(BestIndex,:); 
                 
            currentIteration=currentIteration+1;
        end
        
         BestFitness = min(fitness)
         BestGen

