function D = getElasticMat(analysis, intMode,E,v)

        switch analysis

            case 'PS'
                D{1}=(E/(1-v^2))*[1 v 0; v 1 0; 0 0 (1-v)/2];
            case 'PE'
                if strcmp(intMode,'SRI')
                    Dp=(E/3)*[1/(1-2*v) 1/(1-2*v)      0;
                          1/(1-2*v) 1/(1-2*v)      0;
                          0         0         3/(2*(1+v))];

                    Ds=E/(3*(1+v))*[2 -1 0;
                                   -1  2 0;
                                    0  0 0];
                    D{1}=Dp;
                    D{2}=Ds;
                else
                    D{1}=(E/((1+v)*(1-2*v)))*[1-v v 0; v 1-v 0; 0 0 (1-2*v)/2];
                end              

        end


    end