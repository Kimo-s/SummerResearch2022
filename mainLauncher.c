#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> 
#include <unistd.h> /* for fork */
#include <sys/types.h> /* for pid_t */
#include <sys/wait.h> /* for wait */

int main(int argc, char *argv[]){
    int numOfInsteances = 10;
    if(argc == 2){
        int temp = atoi(argv[1]);
        if(temp == 0) {
            return  0;
        }
        numOfInsteances = temp;
    }

    
    pid_t pids[numOfInsteances];
    /*Spawn a child to run the program.*/

    srand(time(0));
    int seed;
    for(int i = 0; i < numOfInsteances; i++){
        seed = rand();
        pids[i] = fork();

        //Child process
        if(pids[i] == 0){
            char snum[256];
            sprintf(snum, "%d", seed);
            char *argv[]={"sim",snum,NULL};
            execv("sim",argv);
            exit(127); /* only if execv fails */
        }
    }

    for(int i = 0; i < numOfInsteances; i++){
        waitpid(pids[i],0,0);
    }



    // pid_t pid = fork();
    // if (pid==0) { /* child process */
    //     static char *argv[]={"sim","123123122",NULL};
    //     execv("sim",argv);
    //     exit(127); /* only if execv fails */
    // }
    // else { /* pid!=0; parent process */
    //     waitpid(pid,0,0); /* wait for child to exit */
    // }

}