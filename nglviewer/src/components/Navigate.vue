<template>
    <div class="navigation-container">
        <button
            class="nav-button"
            :disabled="currentSample <= 1"
            @click="prevSample"
        >
            <i class="fas fa-arrow-left"></i>
        </button>

        <div class="sample-counter">{{ currentSample }}/{{ totalSamples }}</div>

        <button
            class="nav-button"
            :disabled="currentSample >= totalSamples"
            @click="nextSample"
        >
            <i class="fas fa-arrow-right"></i>
        </button>
        <div v-if="showMinimizeCheckbox" class="minimize-checkbox">
            <input
                type="checkbox"
                id="minimized"
                v-model="minimizeView"
                @change="handleMinimizeChange"
            />
            <label for="minimized">Minimized</label>
        </div>
    </div>
</template>

<script>
import axios from "axios";

export default {
    props: {
        processId: {
            type: String,
            default: null,
        },
        totalSamples: {
            type: Number,
            required: true,
        },
        minimize: {
            type: Number,
            default: 0,
        },
    },

    data() {
        return {
            currentSample: 1,
            minimizeView: false, // track the minimize checkbox
        };
    },

    mounted() {
        // Load initial data
        if (this.processId) {
            this.fetchPdbData(this.currentSample);
        }
    },

    computed: {
        showMinimizeCheckbox() {
            return this.currentSample <= this.minimize;
        },
    },

    methods: {
        async fetchPdbData(sampleNumber, minimized = false) {
            if (!this.processId) return;
            try {
                let url;
                if (minimized && this.minimize != 0) {
                    const paddedNumber = String(sampleNumber - 1).padStart(
                        4,
                        "0",
                    );
                    url = `/api/pdb?Result/${this.processId}/minimize/minimized_${paddedNumber}.pdb`;
                } else {
                    const paddedNumber = String(sampleNumber - 1).padStart(
                        4,
                        "0",
                    );
                    url = `/api/pdb?Result/${this.processId}/sample/sample_${paddedNumber}.pdb`;
                }
                const response = await axios.get(url, {
                    headers: {
                        "Cache-Control": "no-cache",
                        Pragma: "no-cache",
                    },
                });
                if (response.data) {
                    this.$emit("update-sequence", response.data);
                }
            } catch (error) {
                console.error("Error fetching PDB data:", error);
            }
        },
        async prevSample() {
            if (this.currentSample > 1) {
                this.currentSample--;
                await this.fetchPdbData(this.currentSample, this.minimizeView);
            }
        },
        async nextSample() {
            if (this.currentSample < this.totalSamples) {
                this.currentSample++;
                await this.fetchPdbData(this.currentSample, this.minimizeView);
            }
        },
        handleMinimizeChange() {
            this.fetchPdbData(this.currentSample, this.minimizeView);
        },
    },
};
</script>

<style scoped>
.navigation-container {
    position: fixed;
    bottom: 2rem;
    left: 50%;
    transform: translateX(-50%);
    display: flex;
    align-items: center;
    gap: 2rem;
    background-color: var(--background-color);
    padding: 0.5rem 1rem;
    border-radius: 25px;
    box-shadow: 0 2px 6px rgba(0, 0, 0, 0.1);
    z-index: 1000;
}

.nav-button {
    background: none;
    border: none;
    font-size: 1rem;
    color: var(--text-color);
    cursor: pointer;
    padding: 0.3rem;
    transition: opacity 0.2s ease;
}

.nav-button:disabled {
    opacity: 0.3;
    cursor: not-allowed;
}

.nav-button:not(:disabled):hover {
    color: var(--primary-color);
}

.sample-counter {
    font-size: 0.9rem;
    font-weight: bold;
    color: var(--text-color);
    min-width: 60px;
    text-align: center;
}

.minimize-checkbox {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    color: var(--text-color);
    font-size: 0.9rem;
    white-space: nowrap;
}

.minimize-checkbox input[type="checkbox"] {
    cursor: pointer;
}
</style>
