classdef ProgressText
    
    %% Public methods
   
    methods (Access = public)
        
        function this = ProgressText(progress_steps)
            
            this.N = floor(progress_steps);
            
            if (this.N<0), error('The number of progress steps in the status bar must a positive integer'); end
   
            this.n = 0;
            this.progress = 0;
            
            this.step = 100 / this.N;

        end
        
        function init(this)
            
            fprintf('progress: [%s] %4.1f%%', repmat(' ',1,this.N), this.progress);
            
        end
        
        function update(this, prog)
            
            if (prog < 0 || prog > 100), error('The progress must be in the interval [0 100]'); end
            
            this.progress = prog;
            
            this.n = floor(prog/this.step);
            
            fprintf( [repmat('\b',1, this.N+2+5 ) repmat('=',1,this.n) repmat(' ',1,this.N-this.n) '] ' sprintf('%4.1f',prog) '%%'] )
            
%             if (prog >= this.n*this.step)
%                 n_prev = this.n;
%                 this.n = floor(prog/this.step); 
%                 fprintf( [repmat('\b',1, this.N-n_prev+2+5 ) repmat('=',1,this.n-n_prev) repmat(' ',1,this.N-this.n) '] ' sprintf('%4.1f',prog) '%%'] )
%             else
%                 fprintf( [repmat('\b',1,5) sprintf('%4.1f',prog) '%%'] );
%             end
    
        end
        
        function printInNewLine(this)
            
            fprintf(['\nprogress: [' repmat('=',1,this.n) repmat(' ',1,this.N-this.n) '] ' sprintf('%4.1f',this.progress) '%%']);
            
        end
        
        
    end
    
    %% Static public methods
    
    methods (Access = public, Static)
        
        
        function test()
            
            tic
            prog_text = ProgressText(50);

            prog_text.init();

            pause(0.2);

            m = 50;
            for i=1:m

                pause(0.01);

                prog = i/m * 100;

                prog_text.update(prog);
                
                if (i==25), prog_text.printInNewLine(); end


            end
            fprintf('\n');
            
            toc
            
        end

        
        
    end
    
    %% Private properties
    
    properties (Access = private)
    
        N % number of discrete steps in status bar
        n % number of '=', representing the progress
        
        step % numeric value of each '='
        
        progress % numeric value of progress (%)
    
    end
end