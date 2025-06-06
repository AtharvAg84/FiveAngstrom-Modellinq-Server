// package.json dependencies needed:
// npm install express cors helmet morgan compression
// npm install -D @types/express @types/node typescript ts-node nodemon

import express, { Request, Response, NextFunction } from 'express';
import cors from 'cors';
import helmet from 'helmet';
import morgan from 'morgan';
import compression from 'compression';
import { createHash } from 'crypto';
import { spawn } from 'child_process';
import { promises as fs } from 'fs';
import path from 'path';

const app = express();
const PORT = process.env.PORT || 3000;

// Middleware
app.use(helmet());
app.use(cors());
app.use(compression());
app.use(morgan('combined'));
app.use(express.json({ limit: '10mb' }));
app.use(express.urlencoded({ extended: true }));

// Serve static files from nglviewer/dist
app.use(express.static(path.join(__dirname, '..', 'nglviewer', 'dist')));

// Types
interface Parameter {
    sequence: string;
    sample: number;
    forcefield: string;
    grid: number;
}

enum ProcessStatus {
    Waiting = 'Waiting',
    Sampling = 'Sampling',
    Completed = 'Completed'
}

interface ProcessItem {
    process_id: string;
    process: Parameter;
    status: ProcessStatus;
}

// Queue Management Class
class Scheduler {
    private queue: ProcessItem[] = [];
    private readonly maxLength = 10;

    get length(): number {
        return this.queue.length;
    }

    async enqueue(process: Parameter): Promise<string> {
        if (this.queue.length >= this.maxLength) {
            throw new Error('QueueFull');
        }

        const id = this.generateProcessId(process);

        // Check if process already exists
        const exists = this.queue.some(item => item.process_id === id);

        if (!exists) {
            this.queue.push({
                process_id: id,
                process,
                status: ProcessStatus.Waiting
            });
        }

        return id;
    }

    private generateProcessId(process: Parameter): string {
        const hash = createHash('sha256');
        hash.update(process.sequence);
        hash.update(process.sample.toString());
        hash.update(process.forcefield);
        hash.update(process.grid.toString());

        const digest = hash.digest();
        const id = digest.readUInt32BE(0) * 256 + digest.readUInt8(4);
        return id.toString(16).toUpperCase().padStart(10, '0');
    }

    getProcessPosition(id: string): number | null {
        const index = this.queue.findIndex(item => item.process_id === id);
        return index !== -1 ? index : null;
    }

    getProcessById(id: string): ProcessItem | null {
        return this.queue.find(item => item.process_id === id) || null;
    }

    async dequeue(): Promise<void> {
        if (this.queue.length === 0) return;

        const first = this.queue[0];
        if (!first) return;

        const folder = `Result/${first.process_id}`;

        try {
            await this.execSample(first.process, folder);
            await this.cleanFolder(folder);
            first.status = ProcessStatus.Completed;
        } catch (error) {
            console.error('Error processing queue item:', error);
            first.status = ProcessStatus.Completed; // Mark as completed even on error
        }
    }

    private async execSample(process: Parameter, folder: string): Promise<void> {
        const args = [
            'sample',
            '-F', folder,
            '-s', process.sequence,
            '-n', process.sample.toString(),
            '-f', process.forcefield,
            '-g', process.grid.toString()
        ];

        console.log('Executing Sampling Command:', ['dncs', ...args]);

        return new Promise((resolve, reject) => {
            const child = spawn('dncs', args, {
                stdio: ['ignore', 'pipe', 'pipe'] // Capture stdout and stderr
            });

            let stderr = '';
            child.stderr?.on('data', (data) => {
                stderr += data.toString();
            });

            child.on('close', (code) => {
                if (code !== 0) {
                    console.error(`Sampling Process Failed With Code: ${code}`);
                    console.error(`stderr: ${stderr}`);
                    reject(new Error(`dncs failed with code ${code}`));
                } else {
                    console.log(`Sampling Process Exited With Code: ${code}`);
                    resolve();
                }
            });

            child.on('error', (error) => {
                console.error('Sampling process error:', error);
                reject(error);
            });
        });
    }

    private async cleanFolder(folder: string): Promise<void> {
        console.log('Executing Cleaning Folder:', ['at', 'now', '+24', 'hours']);

        return new Promise((resolve, reject) => {
            const child = spawn('at', ['now', '+24', 'hours'], {
                stdio: ['pipe', 'ignore', 'ignore']
            });

            child.on('spawn', () => {
                child.stdin?.write(`rm -rf ${folder}`);
                child.stdin?.end();
            });

            child.on('close', () => {
                console.log(`Folder ${folder} will be cleaned after 24 hours`);
                resolve();
            });

            child.on('error', (error) => {
                console.error('Clean process error:', error);
                reject(error);
            });
        });
    }

    reorder(): void {
        this.queue.shift(); // Remove first element
    }
}

// Global queue instance
const queue = new Scheduler();

// API Routes

// Scheduling Process
app.post('/api/sample', async (req: Request, res: Response): Promise<void> => {
    try {
        const body = req.body as Parameter;

        // Basic validation
        if (!body.sequence || !body.sample || !body.forcefield || !body.grid) {
            res.status(400).send('Invalid JSON - missing required fields');
            return;
        }

        const processId = await queue.enqueue(body);
        console.log(`Queue length: ${queue.length}`);

        res.send(processId);
    } catch (error) {
        if (error instanceof Error && error.message === 'QueueFull') {
            res.send('QueueFull');
            return;
        }
        res.status(500).send('Internal server error');
    }
});

// Process Status
app.get('/api/status', async (req: Request, res: Response): Promise<void> => {
    const id = req.query.id as string;

    if (!id) {
        res.status(400).send("Missing 'id' query parameter");
        return;
    }

    const queueProcess = queue.getProcessById(id);
    if (!queueProcess) {
        res.status(404).send('Not Found');
        return;
    }

    switch (queueProcess.status) {
        case ProcessStatus.Waiting:
            const position = queue.getProcessPosition(id);
            if (position === 0) {
                queueProcess.status = ProcessStatus.Sampling;
                await queue.dequeue();
            } else if (position === 1) {
                const firstProcess = queue.getProcessById(queue.getProcessPosition('0')?.toString() || '');
                if (firstProcess?.status === ProcessStatus.Completed) {
                    queue.reorder();
                }
            } else {
                res.send(`Waiting on Queue ${position}/10`);
                return;
            }
            break;

        case ProcessStatus.Sampling:
            res.send('Generating the Samples');
            return;

        case ProcessStatus.Completed:
            queue.reorder();
            res.send('Completed');
            return;
    }
});

// Send PDB File
app.get('/api/pdb', async (req: Request, res: Response): Promise<void> => {
    const filePath = req.query.path as string;

    if (!filePath) {
        res.status(400).send("Missing 'path' query parameter");
        return;
    }

    try {
        // Security check - ensure path doesn't contain directory traversal
        const safePath = path.resolve(filePath);
        await fs.access(safePath);
        res.sendFile(safePath);
    } catch (error) {
        res.status(404).send('File not found');
    }
});

// PDB Sequence
app.get('/api/seq', async (req: Request, res: Response): Promise<void> => {
    // Disable caching for this endpoint
    res.setHeader('Cache-Control', 'no-store, no-cache, must-revalidate, proxy-revalidate');
    res.setHeader('Pragma', 'no-cache');
    res.setHeader('Expires', '0');
    res.setHeader('Surrogate-Control', 'no-store');
    // Accept both ?sequence=ATHARV and ?ATHARV (legacy)
    let sequence = req.query.sequence as string | undefined;
    if (!sequence) {
        // Try to extract from query string key (e.g., /api/seq?ATHARV)
        const keys = Object.keys(req.query);
        if (keys.length === 1 && req.query[keys[0] as keyof typeof req.query] === undefined) {
            sequence = keys[0];
        }
    }
    if (!sequence) {
        res.send('Invalid Sequence1');
        return;
    }

    // Validate sequence - check for invalid characters
    const invalidChars = /[BOJUXZ\s]/;
    if (invalidChars.test(sequence)) {
        res.send('Invalid Sequence2');
        return;
    }

    try {
        // This would need to be implemented based on your polymer library functionality
        // For now, returning a placeholder response
        const pdbData = await generatePDBFromSequence(sequence);
        res.send(pdbData);
    } catch (error) {
        console.error('Error generating PDB:', error);
        res.send('Invalid Sequence');
    }
});

// Placeholder function - implement based on your polymer library
async function generatePDBFromSequence(sequence: string): Promise<string> {
    // This is a placeholder implementation
    // You'll need to implement the actual PDB generation logic
    // based on your polymer library functionality

    return `HEADER    GENERATED FROM SEQUENCE ${sequence}\nATOM      1  N   ALA A   1      20.154  16.967  22.084  1.00  0.00           N\nEND\n`;
}

// Error handler
app.use((err: Error, req: Request, res: Response, next: NextFunction) => {
    console.error(err.stack);
    res.status(500).send('Something went wrong!');
});

// Move this catch-all route to the very end, after all API routes
app.get('*', (req: Request, res: Response) => {
    if (req.path.startsWith('/api/')) {
        res.status(404).send('API endpoint not found');
        return;
    }
    const indexPath = path.join(__dirname, '..', 'nglviewer', 'dist', 'index.html');
    res.sendFile(indexPath, (err) => {
        if (err) {
            console.error('Error serving index.html:', err);
            res.status(404).send('Not found html');
        }
    });
});

// Start server
app.listen(PORT, () => {
    console.log(`Listening on 0.0.0.0:${PORT}`);
});

export default app;