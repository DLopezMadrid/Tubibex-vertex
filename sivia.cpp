#include "sivia.h"
#include <QList>
#include "vibes.h"
#include <QDateTime>

#define __PREC__ 1e-11
#define __METH__ RK4
#define __DURATION__ 0.5

/// Initializes constraints and prepares data to be processed
Sivia::Sivia(Data &data, bool calcInner)
{
    //constraints to check to see if a trajectory belongs to a tube

    //c1=   gdot= d/dx(gi)(x,t)*f(x,t)+d/dt(gi)(x,t)>=0
    //c2=   gi(x,t)=0j
    //c3=   g(x,t)<=0

    Function g("g.txt");



    Function dg(g, Function::DIFF);                                     //  d/dx(gi)(x,t)
    Variable x(data.numVarF),t;                                         //we have x[] and t as variables for our fns

    // initialize auxMat and auxVector to the correct sizes and fill with zeros
    IntervalMatrix auxMat(data.numVarF+1, data.numVarF,Interval::ZERO);
    IntervalVector auxVector(data.numVarF+1,Interval::ZERO);

    //put 1 in the diagonal of auxMat
    for (int i=0; i<data.numVarF; i++){
        auxMat[i][i]=1;}

    auxVector[data.numVarF]=1;

    if (calcInner){         //for the inner approximation of the tube, we set a new function f and its correspondent constraints
        cout<<endl<<"Start inner approx calculation"<<endl;
        Function f=("f.txt");
        Function gdot(x,t,dg(x,t)*(auxMat*transpose(f(x,t))+auxVector));


        NumConstraint c3(g, EQ);                                           //LEQ means less or equal 0

        Array<Ctc> individualTubeConstraints;                               //put together the constraints in an array

        individualTubeConstraints.add(*new CtcHC4(Array<NumConstraint>(c3)));

        CtcUnion unionTubeConstraints(individualTubeConstraints);           //calculate the union

        Ctc3BCid tubeConstraints(unionTubeConstraints);                     //this contracts the whole union to its final form

        data.boxes.push_back(data.initialBox);                              //initialize the boxes

        do_Sivia(tubeConstraints, data, gdot, calcInner);

        print_results(data);
    }

    else{                       //for the outer approximation
        Function f("f.txt");
        Function gdot(x,t,dg(x,t)*(auxMat*transpose(f(x,t))+auxVector));
        //c1 & c2
        Array<NumConstraint> c1, c2;

        int numConstraints = data.g->expr().dim.max_index()+1;              //find how many gi we have

        for (int i = 0; i < numConstraints; ++i) {                          //create constraints based on the dimensions of g
            c1.add(*new NumConstraint(x,t,gdot(x,t)[i] >= 0));
            c2.add(*new NumConstraint(x,t,g(x,t)[i] = 0));
        }

        NumConstraint c3(g, LEQ);                                           //LEQ means less or equal 0

        Array<Ctc> individualTubeConstraints;                               //put together the constraints in an array

        for (int i=0;i<numConstraints ;i++) {
            individualTubeConstraints.add(*new CtcHC4(Array<NumConstraint>(c1[i],c2[i],c3)));}

        CtcUnion unionTubeConstraints(individualTubeConstraints);           //calculate the union

        Ctc3BCid tubeConstraints(unionTubeConstraints);                     //this contracts the whole union to its final form

        data.boxes.push_back(data.initialBox);                              //initialize the boxes

        do_Sivia(tubeConstraints, data, gdot, calcInner);
        if (!data.calcInner){
            print_results(data);
        }
    }



}



/// Processes the data using contractors and bissections. Classifies the boxes in outside (grey), back_in(yellow) and unsafe (red)
void Sivia::do_Sivia(Ctc& tubeConstraints, Data &data, Function gdot, bool calcInner) {

    QTime tSivia;
    tSivia.start();

    if (calcInner)                  //inner approximation calculation
    {
        int count=0;
        while (!data.boxes.empty()) {
            IntervalVector currentBox = data.boxes.front();                 //start from the first one
            data.boxes.pop_front();                                         //once it has been copied remove the first box

            IntervalVector auxBox=currentBox;                               //store it in aux variable to compare later

            tubeConstraints.contract(currentBox);                           //contract the current box using the previously calculated constraints
            if (currentBox!=auxBox){                                        //if the box has been contracted
                IntervalVector* removedByContractorInner;
                int setDiff=auxBox.diff(currentBox, removedByContractorInner);   //set difference between the contracted box and the original box
                for (int i = 0; i < setDiff; ++i) {
                    //data.boxesOutside.push_back(removedByContractor[i]);    //add the areas removed by the contractor to the outside set


                    bool testInside=true;
                    IntervalVector gg=data.g->eval_vector(removedByContractorInner[i]);

                    for(int j = 0; j<gg.size(); j++){
                        testInside = testInside && (gg[j].ub()<=0);
                    }
                    if (testInside) {
                        data.boxesInside.append(removedByContractorInner[i]);
                    }

                }
                delete[] removedByContractorInner;

            }

            if(data.realTimeDraw){                                          //draw the boxes processing in real time
                draw_update(data, auxBox, currentBox);
            }


            bool allBoxesLessEpsilon=true;                                                                              //check if all the boxes are smaler than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox.diam()[i])<=data.epsilons[i]));
            }
            allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox[currentBox.size()-1].diam())<=data.dt));         //check the time box also


            bool boxesLessEpsilon=false;                                                                                //check if at least one box is smaller than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                boxesLessEpsilon = boxesLessEpsilon||((currentBox[i].diam())<=data.epsilons[i]);
            }
            boxesLessEpsilon = boxesLessEpsilon&&((currentBox[currentBox.size()-1].diam())<=data.dt);                   //check time box



            if (boxesLessEpsilon && !allBoxesLessEpsilon){
                IntervalVector xnext = currentBox.subvector(0, data.numVarF-1).mid();   //using the middle point of the box calculate the future positions using euler method
                IntervalVector x = currentBox.mid();
                bool testBackIn;
                for (int i = 0;i<data.numFuturePos;i++){                                // Euler method: x(n+1)=x(n)+dt*fx
                    x[data.numVarF]= x[data.numVarF].mid();
                    testBackIn = true;
                    xnext=xnext+(data.dt)*data.f->eval_vector(x);
                    x.put(0, xnext);
                    x[data.numVarF] = x[data.numVarF]+(data.dt);
                    IntervalVector gg=data.g->eval_vector(x);
                    for(int j = 0; j<gg.size(); j++){
                        testBackIn = testBackIn && (gg[j].ub()<0);                      //test if it comes back to the bubble

                    }
                    if(testBackIn == true){                                             //If so we calculate the max deviation
                        break;
                    }
                }

                if(testBackIn == true){                                                 //If my box was back in the bubble after integration, I store it in boxesbackin
                    (data.boxesInsideBackIn).append(currentBox);

                    continue;
                }
            }


            if (allBoxesLessEpsilon) {                                                          //if allBoxesLessEpsilon = true the box is unsafe and I continue my loop
                (data.boxesInsideUnsafe).push_back(currentBox);
                count++;
                if (count >=data.maxNumUnsafeBoxes && data.maxNumUnsafeBoxesActivated){         //If I have more boxes than nbPerhaps I stop the loop and I display the results
                    break;
                }
            }
            else {                                                                              //Otherwise we bissect following the widest diameter
                double l = 0;
                double l_temp = 0;
                int v = -1;
                for(int i = 0; i<currentBox.size()-1; i++){                                     //test that the diameter of the boxes doesnt depend on time
                    if(currentBox[i].is_bisectable()||!(currentBox[i].is_degenerated())){
                        l_temp = currentBox[i].diam();
                        if(l_temp>=data.epsilons[i] && l_temp/(data.epsilons[i]) > l){
                            l = l_temp/(data.epsilons[i]);
                            v = i;
                        }
                    }
                }

                l_temp = currentBox[currentBox.size()-1].diam();                                //test the time interval
                if(l_temp>=data.dt && l_temp/(data.dt) > l){
                    v = currentBox.size()-1;
                }
                if(v != -1 && currentBox[v].is_bisectable()){                                   // then the test interval of the state variables, and then it bisects the interval which has the largest diameter
                    pair<IntervalVector,IntervalVector> boxes=currentBox.bisect(v, 0.5);
                    (data.boxes).push_back(boxes.first);
                    (data.boxes).push_back(boxes.second);
                }
                else{
                    if (data.myDebug){
                        std::cout<<"Cannot be bisected \n";
                    }
                }
            }
        }

    }


    else                            //outer approximation
    {
        int count=0;
        //SIVIA
        //process all the boxes in data
        while (!data.boxes.empty()) {
            IntervalVector currentBox = data.boxes.front();                 //start from the first one
            data.boxes.pop_front();                                         //once it has been copied remove the first box

            IntervalVector auxBox=currentBox;                               //store it in aux variable to compare later

            tubeConstraints.contract(currentBox);                           //contract the current box using the previously calculated constraints
            if (currentBox!=auxBox){                                        //if the box has been contracted
                IntervalVector* removedByContractor;
                int setDiff=auxBox.diff(currentBox, removedByContractor);   //set difference between the contracted box and the original box
                for (int i = 0; i < setDiff; ++i) {
                    data.boxesOutside.push_back(removedByContractor[i]);    //add the areas removed by the contractor to the outside set
                }
                delete[] removedByContractor;
            }

            if(data.realTimeDraw){                                          //draw the boxes processing in real time
                draw_update(data, auxBox, currentBox);
            }


            bool allBoxesLessEpsilon=true;                                                                              //check if all the boxes are smaler than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox.diam()[i])<=data.epsilons[i]));
            }
            allBoxesLessEpsilon = (allBoxesLessEpsilon && ((currentBox[currentBox.size()-1].diam())<=data.dt));         //check the time box also


            bool boxesLessEpsilon=false;                                                                                //check if at least one box is smaller than epsilon
            for (int i=0;(i<(currentBox.size()-1));i++){
                boxesLessEpsilon = boxesLessEpsilon||((currentBox[i].diam())<=data.epsilons[i]);
            }
            boxesLessEpsilon = boxesLessEpsilon&&((currentBox[currentBox.size()-1].diam())<=data.dt);                   //check time box


            if (boxesLessEpsilon && !allBoxesLessEpsilon){
                IntervalVector xnext = currentBox.subvector(0, data.numVarF-1).mid();   //using the middle point of the box calculate the future positions using euler method
                IntervalVector x = currentBox.mid();
                bool testBackIn;
                for (int i = 0;i<data.numFuturePos;i++){                                // Euler method: x(n+1)=x(n)+dt*fx
                    x[data.numVarF]= x[data.numVarF].mid();
                    testBackIn = true;
                    xnext=xnext+(data.dt)*data.f->eval_vector(x);
                    x.put(0, xnext);
                    x[data.numVarF] = x[data.numVarF]+(data.dt);
                    IntervalVector gg=data.g->eval_vector(x);
                    for(int j = 0; j<gg.size(); j++){
                        testBackIn = testBackIn && (gg[j].ub()<0);                      //test if it comes back to the bubble

                    }
                    if(testBackIn == true){                                             //If so we calculate the max deviation
                        break;
                    }
                }

                //                if(testBackIn == true){                                                 //If my box was back in the bubble after integration, I store it in boxesbackin
                //                    (data.boxesBackIn).append(currentBox);
                //                    continue;
                //                }
            }


            if (allBoxesLessEpsilon) {                                                          //if allBoxesLessEpsilon = true the box is unsafe and I continue my loop
                (data.boxesUnsafe).push_back(currentBox);
                count++;
                if (count >=data.maxNumUnsafeBoxes && data.maxNumUnsafeBoxesActivated){         //If I have more boxes than nbPerhaps I stop the loop and I display the results
                    break;
                }
            }
            else {                                                                              //Otherwise we bissect following the widest diameter
                double l = 0;
                double l_temp = 0;
                int v = -1;
                for(int i = 0; i<currentBox.size()-1; i++){                                     //test that the diameter of the boxes doesnt depend on time
                    if(currentBox[i].is_bisectable()||!(currentBox[i].is_degenerated())){
                        l_temp = currentBox[i].diam();
                        if(l_temp>=data.epsilons[i] && l_temp/(data.epsilons[i]) > l){
                            l = l_temp/(data.epsilons[i]);
                            v = i;
                        }
                    }
                }

                l_temp = currentBox[currentBox.size()-1].diam();                                //test the time interval
                if(l_temp>=data.dt && l_temp/(data.dt) > l){
                    v = currentBox.size()-1;
                }
                if(v != -1 && currentBox[v].is_bisectable()){                                   // then the test interval of the state variables, and then it bisects the interval which has the largest diameter
                    pair<IntervalVector,IntervalVector> boxes=currentBox.bisect(v, 0.5);
                    (data.boxes).push_back(boxes.first);
                    (data.boxes).push_back(boxes.second);
                }
                else{
                    if (data.myDebug){
                        std::cout<<"Cannot be bisected \n";
                    }
                }
            }
        }




        for(int i=0; i<data.boxesUnsafe.size();i++) {                   //process unsafe boxes

            IntervalVector currentBox=data.boxesUnsafe.at(i);
            IntervalVector nextBox = currentBox.subvector(0, data.numVarF-1);

            if (data.intMethod==0){                                     //Guaranteed integration

                // State variables
                Variable y(data.numVarF);

                // Initial conditions
                IntervalVector yinit(data.numVarF);

                for (int i = 0; i < data.numVarF; ++i) {
                    yinit[i] = currentBox[i];
                    cout<<currentBox[i]<<endl;
                }

                Interval t = currentBox[data.numVarF];
                Interval xd = 7*t;
                Interval xdd = 7;
                Interval yd = sin(0.1*t);
                Interval ydd = 0.1*cos(0.1*t);
                //                Interval xdiff = (xd-y[0]+xdd);
                //                Interval ydiff = (yd-y[1]+ydd);
                //                Interval norm =  (sqrt((xdiff)^2 +(ydiff)^2));


                // system fn has to be re entered here, cannot be loaded directly from text file
                //pendulum
                Function ydot = Function (y,Return (y[1],-1*sin(y[0])-0.15*y[1]));                               // Ivp contruction (initial time is 0.0)



                //non holonomic
                //                Function ydot = Function (y, Return ((sqrt(((xd-y[0]+xdd)*(xd-y[0]+xdd)) +((yd-y[1]+ydd)*(yd-y[1]+ydd))))*cos(y[2]), (sqrt(((xd-y[0]+xdd)*(xd-y[0]+xdd)) +((yd-y[1]+ydd)*(yd-y[1]+ydd))))*sin(y[2]),10*(cos(y[2])*((yd-y[1]+ydd))-sin(y[2])*((xd-y[0]+xdd)))/(sqrt(((xd-y[0]+xdd)*(xd-y[0]+xdd)) +((yd-y[1]+ydd)*(yd-y[1]+ydd))))));                               // Ivp contruction (initial time is 0.0)

                QTime t1;
                t1.start();

                ivp_ode problem = ivp_ode (ydot,0.0 , yinit);

                // Simulation construction and run
                simulation simu = simulation (&problem,data.dt*data.numFuturePos, __METH__, __PREC__);          //uses Runge-kutta4 method
                data.boxesUnsafeFuture.append(simu.run_simulation());                                           //modified ibex_simulation.h to make it return a list with all the solutions, not just the last one

                double timeSiviaCalculations1=t1.elapsed()/1000.0;
                double timeSiviaCalculationsTotal=tSivia.elapsed()/1000.0;
                cout<<endl<<"Unsafe # "<<i<<"  , Box time = "<<timeSiviaCalculations1<<" , Total elapsed time = "<<timeSiviaCalculationsTotal<<endl;
            }

            if (data.intMethod==1){                                     //euler method

                for (int i = 0;i<data.numFuturePos;i++){

                    IntervalVector gdotValue=gdot.eval_vector(currentBox);                                      //evaluate the g and gdot functions to inspect the constraints
                    IntervalVector gValue=data.g->eval_vector(currentBox);

                    if (data.myDebug){
                        cout<<"box = "<<currentBox<<endl;

                        for(int j = 0; j<gValue.size(); j++){
                            cout<<"gdot"<<j<<" = "<<gdotValue<<"  /  "<<(gdotValue[j].lb()>0) <<endl;                       //gdot i values
                        }
                    }


                    nextBox=nextBox+(data.dt)*data.f->eval_vector(currentBox);              //euler method
                    data.boxesUnsafeFuture.append(nextBox);

                    currentBox.put(0, nextBox);
                    currentBox[data.numVarF] = currentBox[data.numVarF]+(data.dt);          //increase time for the next step
                }
            }
            if (data.intMethod==2){                                     //vertex unsafe


                //check derivative of the unsafe boxes of outer g wrt the inner g
                //                Function g_inner("g_inner.txt");

                //                Function dg_inner(g_inner, Function::DIFF);                                     //  d/dx(gi)(x,t)
                //                Variable x(data.numVarF),t;                                         //we have x[] and t as variables for our fns

                //                // initialize auxMat and auxVector to the correct sizes and fill with zeros
                //                IntervalMatrix auxMat(data.numVarF+1, data.numVarF,Interval::ZERO);
                //                IntervalVector auxVector(data.numVarF+1,Interval::ZERO);

                //                //put 1 in the diagonal of auxMat
                //                for (int i=0; i<data.numVarF; i++){
                //                    auxMat[i][i]=1;}

                //                auxVector[data.numVarF]=1;

                //                Function f=("f.txt");
                //                Function gdot_inner(x,t,dg_inner(x,t)*(auxMat*transpose(f(x,t))+auxVector));

                //                cout<<"gdot: "<<gdot_inner<<endl;
                //                IntervalVector gdot_inner_Result=gdot_inner.eval_vector(currentBox);

                //                bool testNegative=true;

                //                for(int j = 0; j<gdot_inner_Result.size(); j++){
                //                    testNegative = testNegative && (gdot_inner_Result[j].ub()<0);
                //                }
                //                if (testNegative) {
                //                    data.boxesOuterGSafeForInnerG.append(currentBox);
                //                    cout<<"Safe for inner G: "<<currentBox<<" result: "<<gdot_inner_Result<<endl;
                //                }
                //                else{
                //                    data.boxesOuterGUnsafeForInnerG.append(currentBox);
                //                    cout<<"Unsafe for inner G: "<<currentBox<<" result: "<<gdot_inner_Result<<endl;
                //                }


                //system
                double x_k_lb[currentBox.size()], x_k_ub[currentBox.size()];
                for (int j = 0; j < currentBox.size(); ++j) {
                    x_k_ub[j]=currentBox[j].ub();
                    x_k_lb[j]=currentBox[j].lb();
                }

                int numVar=currentBox.size()-1;
                int numSamples=data.dt*data.numFuturePos*1000;
                double x_k_uu[numVar][numSamples];
                double x_k_ul[numVar][numSamples];
                double x_k_lu[numVar][numSamples];
                double x_k_ll[numVar][numSamples];


                x_k_uu[0][0]=x_k_ub[0];
                x_k_ll[0][0]=x_k_lb[0];
                x_k_lu[0][0]=x_k_lb[0];
                x_k_ul[0][0]=x_k_ub[0];

                x_k_uu[1][0]=x_k_ub[1];
                x_k_ll[1][0]=x_k_lb[1];
                x_k_lu[1][0]=x_k_ub[1];
                x_k_ul[1][0]=x_k_lb[1];

                double g_outer_uu, g_outer_lu, g_outer_ul, g_outer_ll, g_inner_uu, g_inner_lu, g_inner_ul, g_inner_ll;
                bool unsafeBox=false;
                bool safeBox=false;


                //pendulum system
                for (int j = 1; j < data.dt*data.numFuturePos*1000; ++j) {
                    x_k_uu[1][j]=data.dt/10.0*(-sin(x_k_uu[0][j-1])-0.15*x_k_uu[1][j-1])+x_k_uu[1][j-1];
                    x_k_ll[1][j]=data.dt/10.0*(-sin(x_k_ll[0][j-1])-0.15*x_k_ll[1][j-1])+x_k_ll[1][j-1];
                    x_k_lu[1][j]=data.dt/10.0*(-sin(x_k_lu[0][j-1])-0.15*x_k_lu[1][j-1])+x_k_lu[1][j-1];
                    x_k_ul[1][j]=data.dt/10.0*(-sin(x_k_ul[0][j-1])-0.15*x_k_ul[1][j-1])+x_k_ul[1][j-1];

                    x_k_uu[0][j]=(x_k_uu[1][j-1])*(data.dt/10.0)+x_k_uu[0][j-1];
                    x_k_ll[0][j]=(x_k_ll[1][j-1])*(data.dt/10.0)+x_k_uu[0][j-1];
                    x_k_lu[0][j]=(x_k_lu[1][j-1])*(data.dt/10.0)+x_k_uu[0][j-1];
                    x_k_ul[0][j]=(x_k_ul[1][j-1])*(data.dt/10.0)+x_k_uu[0][j-1];

                    //add discretized sys

                    //check if they leave the outer g
                    g_outer_uu= (x_k_uu[0][j]*x_k_uu[0][j])+(x_k_uu[1][j]*x_k_uu[1][j])-1;
                    g_outer_lu= (x_k_lu[0][j]*x_k_lu[0][j])+(x_k_lu[1][j]*x_k_lu[1][j])-1;
                    g_outer_ul= (x_k_ul[0][j]*x_k_ul[0][j])+(x_k_ul[1][j]*x_k_ul[1][j])-1;
                    g_outer_ll= (x_k_ll[0][j]*x_k_ll[0][j])+(x_k_ll[1][j]*x_k_ll[1][j])-1;

                    //check if the trajectories reenter the inner g
                    if (data.ellipseInner) {
                        g_inner_uu= (x_k_uu[0][j]*x_k_uu[0][j])/0.81+(x_k_uu[1][j]*x_k_uu[1][j])/(0.4*0.4)-1;
                        g_inner_lu= (x_k_lu[0][j]*x_k_lu[0][j])/0.81+(x_k_lu[1][j]*x_k_lu[1][j])/(0.4*0.4)-1;
                        g_inner_ul= (x_k_ul[0][j]*x_k_ul[0][j])/0.81+(x_k_ul[1][j]*x_k_ul[1][j])/(0.4*0.4)-1;
                        g_inner_ll= (x_k_ll[0][j]*x_k_ll[0][j])/0.81+(x_k_ll[1][j]*x_k_ll[1][j])/(0.4*0.4)-1;

                    }
                    else{
                        g_inner_uu= (x_k_uu[0][j]*x_k_uu[0][j])/(0.99*0.99)+(x_k_uu[1][j]*x_k_uu[1][j])/(0.96*0.96)-1;
                        g_inner_lu= (x_k_lu[0][j]*x_k_lu[0][j])/(0.99*0.99)+(x_k_lu[1][j]*x_k_lu[1][j])/(0.96*0.96)-1;
                        g_inner_ul= (x_k_ul[0][j]*x_k_ul[0][j])/(0.99*0.99)+(x_k_ul[1][j]*x_k_ul[1][j])/(0.96*0.96)-1;
                        g_inner_ll= (x_k_ll[0][j]*x_k_ll[0][j])/(0.99*0.99)+(x_k_ll[1][j]*x_k_ll[1][j])/(0.96*0.96)-1;

                    }


                    double myFutureBox[j][2][2];


                    double myMaxPos= qMax(qMax(qMax(x_k_uu[0][j],x_k_ll[0][j]),x_k_lu[0][j]),x_k_ul[0][j]);
                    double myMinPos= qMin(qMin(qMin(x_k_uu[0][j],x_k_ll[0][j]),x_k_lu[0][j]),x_k_ul[0][j]);

                    double myMaxVel= qMax(qMax(qMax(x_k_uu[1][j],x_k_ll[1][j]),x_k_lu[1][j]),x_k_ul[1][j]);
                    double myMinVel= qMin(qMin(qMin(x_k_uu[1][j],x_k_ll[1][j]),x_k_lu[1][j]),x_k_ul[1][j]);

                    myFutureBox[j][0][0]=myMinPos;
                    myFutureBox[j][0][1]=myMaxPos;

                    myFutureBox[j][1][0]=myMinVel;
                    myFutureBox[j][1][1]=myMaxVel;

                    IntervalVector BoxFuture(data.numVarG, myFutureBox[j]);
                    data.boxesVertexFuture.append(BoxFuture);

                    if (true &&(g_outer_uu>=0 || g_outer_lu>=0 || g_outer_ul>=0 || g_outer_ll>=0)){
                        unsafeBox=true;
                        data.boxesUnsafeOuterG.append(currentBox);
                        double myUnsafeFutureBox[j][2][2];
                        for (int k = 0; k < j; ++k) {

                            myMaxPos= qMax(qMax(qMax(x_k_uu[0][j],x_k_ll[0][j]),x_k_lu[0][j]),x_k_ul[0][j]);
                            myMinPos= qMin(qMin(qMin(x_k_uu[0][j],x_k_ll[0][j]),x_k_lu[0][j]),x_k_ul[0][j]);

                            myMaxVel= qMax(qMax(qMax(x_k_uu[1][j],x_k_ll[1][j]),x_k_lu[1][j]),x_k_ul[1][j]);
                            myMinVel= qMin(qMin(qMin(x_k_uu[1][j],x_k_ll[1][j]),x_k_lu[1][j]),x_k_ul[1][j]);

                            myUnsafeFutureBox[j][0][0]=myMinPos;
                            myUnsafeFutureBox[j][0][1]=myMaxPos;

                            myUnsafeFutureBox[j][1][0]=myMinVel;
                            myUnsafeFutureBox[j][1][1]=myMaxVel;

                            IntervalVector BoxUnsafeFuture(data.numVarG, myUnsafeFutureBox[k]);
                            data.boxesUnsafeOuterGFuture.append(BoxUnsafeFuture);
                            //                            cout<<"Unsafe vertex box "<<BoxUnsafeFuture<<endl;
                        }
                        break;
                    }

                    if (j>20 && g_inner_uu<0 && g_inner_lu<0 && g_inner_ul<0 && g_inner_ll<0){
                        safeBox=true;
                        //data.boxesUnsafeOuterG.append(currentBox);
                        double mySafeFutureBox[j][2][2];

                        for (int k = 0; k < j; ++k) {

                            myMaxPos= qMax(qMax(qMax(x_k_uu[0][j],x_k_ll[0][j]),x_k_lu[0][j]),x_k_ul[0][j]);
                            myMinPos= qMin(qMin(qMin(x_k_uu[0][j],x_k_ll[0][j]),x_k_lu[0][j]),x_k_ul[0][j]);

                            myMaxVel= qMax(qMax(qMax(x_k_uu[1][j],x_k_ll[1][j]),x_k_lu[1][j]),x_k_ul[1][j]);
                            myMinVel= qMin(qMin(qMin(x_k_uu[1][j],x_k_ll[1][j]),x_k_lu[1][j]),x_k_ul[1][j]);

                            mySafeFutureBox[j][0][0]=myMinPos;
                            mySafeFutureBox[j][0][1]=myMaxPos;

                            mySafeFutureBox[j][1][0]=myMinVel;
                            mySafeFutureBox[j][1][1]=myMaxVel;

                            IntervalVector BoxSafeFuture(data.numVarG, mySafeFutureBox[k]);
                            data.boxesGuaranteedIntegrationUnsafe.append(BoxSafeFuture);
                            //                            cout<<"Safe vertex box "<<BoxSafeFuture<<endl;
                        }
                        break;
                    }
                }
                if (unsafeBox==false && safeBox==true){
                    data.boxesUncertain.append(currentBox);
                }

            }
        }



        for (int i = 0; i < data.boxesUncertain.size(); ++i) {                              //process uncertain boxes

            IntervalVector uncertainBox=data.boxesUncertain.at(i);
            Variable y(data.numVarF);

            // Initial conditions
            IntervalVector yinit(data.numVarF);

            for (int j = 0; j < data.numVarF; ++j) {
                yinit[j] = uncertainBox[j];
                cout<<uncertainBox[j]<<endl;
            }

            Function ydot = Function (y,Return (y[1],-1*sin(y[0])-0.15*y[1]));                               // Ivp contruction (initial time is 0.0)

            QTime t1;
            t1.start();

            ivp_ode problem = ivp_ode (ydot, 0.0 , yinit);

            // Simulation construction and run
            int prevSize=data.boxesUncertainFuture.size();
            simulation simu = simulation (&problem,data.dt*data.numFuturePos, __METH__, __PREC__);          //uses Runge-kutta4 method
            data.boxesUncertainFuture.append(simu.run_simulation());                                           //modified ibex_simulation.h to make it return a list with all the solutions, not just the last one

            int afterSize=data.boxesUncertainFuture.size();
            for (int k = 0; k < (afterSize-prevSize); ++k) {
                data.boxesUncertainFutureIndex.append(i);
            }
            double timeSiviaCalculations1=t1.elapsed()/1000.0;
            double timeSiviaCalculationsTotal=tSivia.elapsed()/1000.0;
            cout<<endl<<"Unsafe # "<<i<<"  , Box time = "<<timeSiviaCalculations1<<" , Total elapsed time = "<<timeSiviaCalculationsTotal<<endl;

        }
        for (int j = 0; j < data.boxesUncertainFuture.size(); ++j) {
            IntervalVector myBox=data.boxesUncertainFuture.at(j);
            bool testInsideOuterG = true;
            Function g_outer("g_outer.txt");
            IntervalVector gg=g_outer.eval_vector(myBox);
            for(int j = 0; j<gg.size(); j++){
                testInsideOuterG = (testInsideOuterG && (gg[j].ub()<0));
            }
            if (testInsideOuterG==false){
                data.boxesUnsafeOuterGFuture.append(myBox);
                data.boxesUnsafeOuterG.append(data.boxesUncertain.at(data.boxesUncertainFutureIndex.at(j)));
            }
            else{
                data.boxesSafeOuterGFuture.append(myBox);
            }
        }
    }
}

/// Uses Vibes to represent the results
void Sivia::print_results(Data data){
    cout<<"Boxes outer approximation = "<<data.boxesOutside.size()+data.boxesUnsafe.size()<<endl;
    cout<<"Outside boxes = "<<data.boxesOutside.size()<<endl;
    cout<<"Unsafe boxes = "<<data.boxesUnsafe.size()<<endl;
    cout<<"Boxes inner approximation = "<<data.boxesInside.size()+data.boxesInsideUnsafe.size()+data.boxesInsideBackIn.size()<<endl;
    cout<<"Inside boxes = "<<data.boxesInside.size()<<endl;
    cout<<"Inside unsafe boxes = "<<data.boxesInsideUnsafe.size()<<endl;

    //use vibes to show the results
    //    bool showFuturePos=false;
    //    init_scene(data, showFuturePos);
    //    draw_update(data, showFuturePos);
    //    draw_vector_field(data);
    //    vibes::endDrawing();

    //    showFuturePos=true;
    //    init_scene(data, showFuturePos);
    //    draw_update(data, showFuturePos);
    //    draw_vector_field(data);
    //    vibes::endDrawing();

    if (data.saveAllPossibleFigures){
        export_all_images(data);}

}

