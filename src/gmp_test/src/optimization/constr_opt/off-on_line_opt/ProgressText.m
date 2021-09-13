classdef ProgressText
    
    %% Public methods
   
    methods (Access = public)
        
        function this = ProgressText(progress_steps)
        
            this.N = progress_steps;
            this.n = 0;
            
            this.step = 100 / this.N;

        end
        
        function init(this)
            
            fprintf('progress: [%s] %4.1f%%', repmat(' ',1,this.N), 0);
            
        end
        
        function update(this, prog)
            
            if (prog >= this.n*this.step)
                n_prev = this.n;
                this.n = floor(prog/this.step); % n + 1; 
                fprintf( [repmat('\b',1, this.N-n_prev+2+5 ) repmat('=',1,this.n-n_prev) repmat(' ',1,this.N-this.n) '] ' sprintf('%4.1f',prog) '%%'] )
            else
                fprintf( [repmat('\b',1,5) sprintf('%4.1f',prog) '%%'] );
            end
    
        end
        
        
    end
    
    %% Static public methods
    
    methods (Access = public, Static)
        
        
        function test()
            
            prog_text = ProgressText(50);

            prog_text.init();

            pause(0.2);

            m = 50;
            for i=1:m

                pause(0.1);

                prog = i/m * 100;

                prog_text.update(prog);


            end
            fprintf('\n');
            
        end

        
        
    end
    
    %% Private properties
    
    properties (Access = private)
    
        N
        n
        
        step
    
    end
end