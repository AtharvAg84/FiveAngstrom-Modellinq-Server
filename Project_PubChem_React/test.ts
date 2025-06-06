import { exec } from 'child_process';

function runDncs(): Promise<string> {
    return new Promise((resolve, reject) => {
        exec('dncs', (error, stdout, stderr) => {
            if (error) {
                reject(stderr || error.message);
            } else {
                resolve(stdout);
            }
        });
    });
}

// Example usage:
runDncs()
    .then(output => console.log('dncs output:', output))
    .catch(err => console.error('Error running dncs:', err));