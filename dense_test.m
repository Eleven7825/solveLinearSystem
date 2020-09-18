% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

% dense test

% Set parameters if such parameters don't exist in work space
disp('Test for Conjugate Gradient and Conjugate Gradient with Preconditioning')
if ~exist('strt','var');     strt = 3;                        end
if ~exist('stop','var');         stop  = 9;                   end
if ~exist('m','var');        m= 2000;                       end
if ~exist('e','var');        e = 0.01;                        end
if ~exist('showplot','var');   showplot=0;                    end
if ~exist('foldername','var');   foldername='Output';         end
if ~exist('fileformat','var'); fileformat='epsc' ;            end
if ~exist('filename1','var');   filename1= 'relation_N_iterations_dense.eps';end
if ~exist('filename2','var');   filename2= 'relation_condN_iterations_dense.eps';end
if ~exist('repeat','var');   repeat= 3;                       end
fprintf('Will run %d rounds of random square matrices with sizes from %d to %d\n',repeat,2^strt,2^stop);
fprintf('Opening file "Output/Report_dense.txt"... ');
fid = fopen([foldername,'/Report_dense.txt'],'w');
disp('done!')

more off
report_1 = sprintf('Report for test for Conjugate Gradient and Conjugate Gradient with Preconditioning\n\n');
report_1 = [report_1, sprintf('This test generates random dense matrices and \nperforms CG and PCG on them and see how many iterations involves for each matrix.\n\n')];
report_1 = [report_1, sprintf('=========================\n')];

% display the parameters
report_1 = [report_1,sprintf('Parameters:\n\n')];
report_1 = [report_1,sprintf('strt = %d\n',strt)];
report_1 = [report_1,sprintf('stop = %d\n',stop)];
report_1 = [report_1,sprintf('e: tolerance error in ||Ax-b|| < e * ||b||. \n')];
report_1 = [report_1,sprintf('m = %d\n',m)];
report_1 = [report_1,sprintf('e = %f\n',e)];
report_1 = [report_1,sprintf('showplot = %s\n',num2str(showplot))];
report_1 = [report_1,sprintf('foldername = %s\n',foldername)];
report_1 = [report_1,sprintf('filename1 = %s\n',filename1)];
report_1 = [report_1,sprintf('filename2 = %s\n',filename2)];
report_1 = [report_1,sprintf('fileformat = %s\n',fileformat)];
report_1 = [report_1,sprintf('repeat = %d\n\n',repeat)];
report_1 = [report_1,sprintf('Parameters can be changed as variables in the workspace.\n')];
report_1 = [report_1,sprintf('=========================\n\n')];
report_1 = [report_1,sprintf("We have generated ramdom dense matrices with sizes \nfrom 2^%d to 2^%d, and have done the CG and PCG for them respectively.\n\n",strt,stop)];
report_1 = [report_1,sprintf("We have generated ramdom dense matrices by letting \na = randi(n,n) and A = a'*a.\n\n")];
report_1 = [report_1,sprintf('We have done the above process for %d rounds, \nyou can change variable named repeat in workspace for rounds executed.\n\n', repeat)];
report_1 = [report_1,sprintf('e: tolerance error in ||Ax-b|| < e * ||b||. \n')];
report_1 = [report_1,sprintf('m: maximum number of iterations allowed\n')];
report_1 = [report_1,sprintf('A number of %d iterations means that the maximum was \nreached and the method was aborted\n\n',m)];
report_1 = [report_1,sprintf('showplot: 1 to show plot window, 0 to hide\n\n')];
report_1 = [report_1,sprintf('Below is iterations CG and P_CG with ramdom matrices\n with sizes N from 2^%d to 2^%d\n',strt,stop)];


% if the folder doesn't exist, create the folder
if ~exist(foldername,'dir'); mkdir(foldername); end

% Get the result from ploting
N = 2.^(strt:stop);
N_1 = [];
x1 = [];
x2 = [];
condNs = [];
for j = 1:repeat
    fprintf('Round %d/%d',j,repeat)
    report_1 = [report_1,sprintf('\n %d rounds / %d rounds in total:\n',j,repeat)];
    report_1 = [report_1,sprintf('cond                    size              CG round          P_CG round\n')];
    for i = N
        a = randi(i,i);
        A = a' * a;
        x0 = rand(i,1);
        b = rand(i,1);
        [~,i2] = CG(A,b,x0,e,m);
        [~,i3] = P_CG(A,b,x0,e,m);
        
        condN = cond(A);

        if isinf(condN)
            fprintf('%6d',i)
            report_1 = [report_1,sprintf('A is not positive definite!\n')];
            continue
        end
        
        N_1 = [N_1,i];
        condNs = [condNs,condN];
        x1 = [x1,i2];
        x2 = [x2,i3];
        report_1 = [report_1,sprintf('%9d            %4d               %5d               %5d\n',condN,i,i2,i3)];
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
legend({'Conjugate Gradient','Preconditioned Conjugate Gradient'},'location','northwest');
disp('done!')

% save the file if the user choose to save it
fprintf(['Saving plot as "' foldername '/' filename1 '"... '])
saveas(fig1,[foldername '/' filename1 ],fileformat);
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

% plot the relationship between condition numbers and iteraions
loglog(condNs,x1,'.',condNs,x2,'.')
title(sprintf('Rounds for iteration until error is less than %f',e))
xlabel('Condition numbers')
ylabel('iteration rounds')
legend({'Conjugate Gradient','Preconditioned Conjugate Gradient'},'location','northwest');
disp('done!')

fprintf(['Saving plot as "' foldername '/' filename2 '"... '])
saveas(fig2,[foldername '/' filename2],fileformat);
fprintf('done!\n')

report_1 = [report_1,sprintf('\n\nThe program finished running %d rounds of tests for CG and PCG on dense matricess.\nbye!\n\n',repeat)];
for text = report_1
    fprintf(fid,text);
end

fprintf('Closing file "Output/Report_dense.txt"... ')
fclose(fid);
fprintf('done!\n')
disp('Bye!')