% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

% runtest

% Set parameters if such parameters don't exist in work space
disp('Test for Steepest Descent and Conjugate Gradient')
if ~exist('strt','var');     strt = 3;                        end
if ~exist('stop','var');         stop  = 6;                   end
if ~exist('m','var');        m= 400000;                       end
if ~exist('e','var');        e = 0.01;                        end
if ~exist('showplot','var');   showplot=0;                    end
if ~exist('foldername','var');   foldername='Output';         end
if ~exist('fileformat','var'); fileformat='epsc' ;            end
if ~exist('filename5','var');   filename5= 'relation_N_iterations_runtest.eps';end
if ~exist('filename6','var');   filename6= 'relation_condN_iterations_dense_runtest.eps';end
if ~exist('repeat','var');   repeat= 3;                       end
fprintf('Will run %d rounds of random square matrices with sizes from %d to %d\n',repeat,2^strt,2^stop);
fprintf('Opening file "Output/Report.txt"... ');

% if the folder doesn't exist, create the folder
if ~exist(foldername,'dir'); mkdir(foldername); end
fid = fopen([foldername,'/Report.txt'],'w');
disp('done!')

more off
report_0 = sprintf('Report for test for Conjugate Gradient and Conjugate Gradient with Preconditioning\n\n');
report_0 = [report_0, sprintf('This test generates random dense matrices and \nperforms CG and PCG on them and see how many iterations involves for each matrix.\n\n')];
report_0 = [report_0, sprintf('=========================\n')];

% display the parameters
report_0 = [report_0,sprintf('Parameters:\n\n')];
report_0 = [report_0,sprintf('strt = %d\n',strt)];
report_0 = [report_0,sprintf('stop = %d\n',stop)];
report_0 = [report_0,sprintf('e: tolerance error in ||Ax-b|| < e * ||b||. \n')];
report_0 = [report_0,sprintf('m = %d\n',m)];
report_0 = [report_0,sprintf('e = %f\n',e)];
report_0 = [report_0,sprintf('showplot = %s\n',num2str(showplot))];
report_0 = [report_0,sprintf('foldername = %s\n',foldername)];
report_0 = [report_0,sprintf('filename5 = %s\n',filename5)];
report_0 = [report_0,sprintf('filename6 = %s\n',filename6)];
report_0 = [report_0,sprintf('fileformat = %s\n',fileformat)];
report_0 = [report_0,sprintf('repeat = %d\n\n',repeat)];
report_0 = [report_0,sprintf('Parameters can be changed as variables in the workspace.\n')];
report_0 = [report_0,sprintf('=========================\n\n')];
report_0 = [report_0,sprintf("We have generated ramdom dense matrices with sizes \nfrom 2^%d to 2^%d, and have done the SD and CG for them respectively.\n\n",strt,stop)];
report_0 = [report_0,sprintf("We have generated ramdom dense matrices by letting \na = randi(n,n) and A = a'*a.\n\n")];
report_0 = [report_0,sprintf('We have done the above process for %d rounds, \nyou can change variable named repeat in workspace for rounds executed.\n\n', repeat)];
report_0 = [report_0,sprintf('e: tolerance error in ||Ax-b|| < e * ||b||. \n')];
report_0 = [report_0,sprintf('m: maximum number of iterations allowed\n')];
report_0 = [report_0,sprintf('A number of %d iterations means that the maximum was \nreached and the method was aborted\n\n',m)];
report_0 = [report_0,sprintf('showplot: 1 to show plot window, 0 to hide\n\n')];
report_0 = [report_0,sprintf('Below is iterations CG and P_CG with ramdom matrices\n with sizes N from 2^%d to 2^%d\n',strt,stop)];


% Get the result from ploting
N = 2.^(strt:stop);
N_1 = [];
x1 = [];
x2 = [];
condNs = [];
for j = 1:repeat
    fprintf('Round %d/%d',j,repeat)
    report_0 = [report_0,sprintf('\n %d rounds / %d rounds in total:\n',j,repeat)];
    report_0 = [report_0,sprintf('cond                    size              SD round          CG round\n')];
    for i = N
        a = randi(i,i);
        A = a' * a;
        x0 = rand(i,1);
        b = rand(i,1);
        [~,i2] = SD(A,b,x0,e,m);
        [~,i3] = CG(A,b,x0,e,m);
        
        condN = cond(A);

        if isinf(condN)
            disp('condition number is Inf!')
            continue
        end
        
        if i2 == m || i3 == m
            fprintf('%6d',i)
            report_0 = [report_0,sprintf('%9d            %4d               %6d               %5d\n',condN,i,i2,i3)];
            continue
        end
        
        N_1 = [N_1,i];
        condNs = [condNs,condN];
        x1 = [x1,i2];
        x2 = [x2,i3];
        report_0 = [report_0,sprintf('%9d            %4d               %6d               %5d\n',condN,i,i2,i3)];
        fprintf('%6d',i)
    end
     fprintf('\n')
end
fprintf('Finish tests.\n\n')

% set up for the ploting
fprintf('Begin plot size N vs the iterations for CG and SD...')
if showplot == 0
    fig1 = figure('visible','off');
elseif showplot == 1
    fig1 = figure('visible','on');
else
    fig1 = figure();
    fprintf('Warning: variable showplot must be 0 or 1.\n')
end

% plot the relationship between N and iteraions
loglog(N_1,x1,'.',N_1,x2,'.')
title(sprintf('Rounds for iteration until error is less than %f',e))
xlabel('size of the matrix')
ylabel('iteration rounds')
legend({'Steepest Descent','Conjugate Gradient'})
disp('done!')

% save the file if the user choose to save it
fprintf(['Saving plot as "' foldername '/' filename5 '"... '])
saveas(fig1,[foldername '/' filename5 ],fileformat);
disp('done!')

% set up for the ploting
fprintf('Begin plot condition numbers vs the iterations for CG and SD...')
if showplot == 0
    fig2 = figure('visible','off');
elseif showplot == 1
    fig2 = figure('visible','on');
else
    fig2 = figure();
    fprintf('Warning: variable showplot must be 0 or 1.\n')
end

% for 0 in x2, force them equal to 1
i = 1;
for x = x2
    if x == 0
        x2(i)=1;
    end
    i = i+1;
end

% theoretical boudries
K = logspace(1,stop+1);
sqK = sqrt(K);
I = log(e/2)./(log((sqK-1)./(sqK+1)));
I2 = log(e)./log((K-1)./(K+1));


% plot the relationship between condition numbers and iteraions
loglog(condNs,x1,'.',condNs,x2,'.',K,I,'-',K,I2,'-')
title(sprintf('Rounds for iteration until error is less than %f',e))
xlabel('Condition numbers')
ylabel('iteration rounds')
legend({'Steepest Descent','Conjugate Gradient'})
disp('done!')

fprintf(['Saving plot as "' foldername '/' filename6 '"... '])
saveas(fig2,[foldername '/' filename6],fileformat);
fprintf('done!\n')

report_0 = [report_0,sprintf('\n\nThe program finished running %d rounds of tests for SDs and CG.\nbye!\n\n',repeat)];
for text = report_0
    fprintf(fid,text);
end

fprintf('Closing file "Output/Report.txt"... ')
fclose(fid);
fprintf('done!\n')
disp('Bye!')