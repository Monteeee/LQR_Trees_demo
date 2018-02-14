# LQR_Trees_demo

Going to be some demo of LQR-Trees control here.
(  Based on the paper:
Tedrake, Russ, et al. "LQR-trees: Feedback motion planning via sums-of-squares verification." The International Journal of Robotics Research 29.8 (2010): 1038-1052. )

System Requirements:

1. Robotics System Toolbox of Matlab2017b.

2. Control System Toolbox

3. SOSTOOLS, release v3.01 (1st July 2016), http://www.cds.caltech.edu/sostools/

<s>4. For SOSTOOLS, the SDPT3 toolbox is needed. The link of compatible SDPT3 toolbox is listed in "System Requirements" of SOSTOOLS page.</s>

4. For SOSTOOLS, install the sedumi toolbox as solver. The link of it is listed in "System Requirements" of SOSTOOLS page. Just download the latest version, extract and copy the folder to matlab toolbox directory.

**BUG fixing for SOSTOOLS:**
In SOSTOOLS, when creating new SOS variable,the automatically generated coefficient has a dimension different from the monomials. This cause the .* operation between them to be ilegal.
To fix this, a feasible but maybe not safe way is to change "V = (coeff.') * ZZSym;" in sosvar, (around line 108 - 110) to:

        if(length(coeff) ~= length(ZZSym))
            temporary = reshape(ZSym * ZSym.', length(ZSym)^2, 1);
            V = (coeff.') * temporary;
        else
            V = (coeff.') * ZZSym;
        end
