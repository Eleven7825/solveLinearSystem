% Copyright (c) 2020 Shiyi Chen and Leonardo T. Rolla
% You can use, modify and redistribute this program under the terms of 
% the GNU Lesser General Public License, either version 3 of the License, 
% or any later version.

% sparse test

% Set parameters if such parameters don't exist in work space
disp('Test for Conjugate Gradient and Conjugate Gradient with Preconditioning')
if ~exist('strt','var');     strt = 3;                        end
if ~exist('stop','var');         stop  = 11;                   end
if ~exist('m','var');        m= 2000;                       end
if ~exist('e','var');        e = 0.01;                        end
if ~exist('showplot','var');   showplot=0;                    end
if ~exist('foldername','var');   foldername='Output';         end
if ~exist('fileformat','var'); fileformat='epsc' ;            end
if ~exist('filename3','var');   filename3= 'relation_N_iterations_sparse.eps';end
if ~exist('filename4','var');   filename4= 'relation_condN_iterations_sparse.eps';end
if ~exist('repeat','var');   repeat= 3;                       end
fprintf('Will run %d rounds of random square matrices with sizes from %d to %d\n',repeat,2^strt,2^stop);
fprintf('Opening file "Output/Report_sparse.txt"... ');
fid = fopen([foldername,'/Report_sparse.txt'],'w');
disp('done!')

more off
report_2 = sprintf('Report for test for Conjugate Gradient and Conjugate Gradient with Preconditioning\n\n');
report_2 = [report_2, sprintf('This test generates random sparse matrices and \nperforms CG and PCG on them and see how many iterations involves for each matrix.\n\n')];
report_2 = [report_2, sprintf('=========================\n')];

% display the parameters
report_2 = [report_2,sprintf('Parameters:\n\n')];
report_2 = [report_2,sprintf('strt = %d\n',strt)];
report_2 = [report_2,sprintf('stop = %d\n',stop)];
report_2 = [report_2,sprintf('e: tolerance error in ||Ax-b|| < e * ||b||. \n')];
report_2 = [report_2,sprintf('m = %d\n',m)];
report_2 = [report_2,sprintf('e = %f\n',e)];
report_2 = [report_2,sprintf('showplot = %s\n',num2str(showplot))];
report_2 = [report_2,sprintf('foldername = %s\n',foldername)];
report_2 = [report_2,sprintf('filename3 = %s\n',filename3)];
report_2 = [report_2,sprintf('filename4 = %s\n',filename4)];
report_2 = [report_2,sprintf('fileformat = %s\n',fileformat)];
report_2 = [report_2,sprintf('repeat = %d\n\n',repeat)];
report_2 = [report_2,sprintf('Parameters can be changed as variables in the workspace.\n')];
report_2 = [report_2,sprintf('=========================\n\n')];
report_2 = [report_2,sprintf("We have generated ramdom sparse matrices with sizes \nfrom 2^%d to 2^%d, and have done the CG and PCG for them respectively.\n\n",strt,stop)];
report_2 = [report_2,sprintf("We have generated ramdom sparse matrices by letting \na = randi(n,n) and A = a'*a, then let all entries be zero except ones \non diagonal and antidiagonal.\n\n")];
report_2 = [report_2,sprintf('We have done the above process for %d rounds, \nyou can change variable named repeat in workspace for rounds executed.\n\n', repeat)];
report_2 = [report_2,sprintf('e: tolerance error in ||Ax-b|| < e * ||b||. \n')];
report_2 = [report_2,sprintf('m: maximum number of iterations allowed\n')];
report_2 = [report_2,sprintf('A number of %d iterations means that the maximum was \nreached and the method was aborted\n\n',m)];
report_2 = [report_2,sprintf('showplot: 1 to show plot window, 0 to hide\n\n')];
report_2 = [report_2,sprintf('Below is iterations CG and P_CG with ramdom matrices\n with sizes N from 2^%d to 2^%d\n',strt,stop)];


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
    report_2 = [report_2,sprintf('\n %d rounds / %d rounds in total:\n',j,repeat)];
    report_2 = [report_2,sprintf('cond                    size              CG round          P_CG round\n')];
    for i = N
        a = randi(i,i);
        A = a' * a;
        for nrow = 1 : i
            for ncol = 1 : i
                if ncol == nrow || (ncol+nrow == i)
                    continue
                else
                    A(nrow,ncol) = 0;
                end
            end
        end
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
        report_2 = [report_2,sprintf('%9d            %4d               %3d               %3d\n',condN,i,i2,i3)];
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
legend({'Conjugate Gradient','Preconditioned Conjugate Gradient'},'location','southeast')
disp('done!')

% save the file if the user choose to save it
fprintf(['Saving plot as "' foldername '/' filename3 '"... '])
saveas(fig1,[foldername '/' filename3 ],fileformat);
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
legend({'Conjugate Gradient','Preconditioned Conjugate Gradient'})
disp('done!')

fprintf(['Saving plot as "' foldername '/' filename4 '"... '])
saveas(fig2,[foldername '/' filename4],fileformat);
fprintf('done!\n')

report_2 = [report_2,sprintf('\n\nThe program finished running %d rounds of tests for CG and PCG on sparse matricess.\nbye!\n\n',repeat)];
for text = report_2
    fprintf(fid,text);
end

fprintf('Closing file "Output/Report_sparse.txt"... ')
fclose(fid);
fprintf('done!\n')

disp('Bye!')