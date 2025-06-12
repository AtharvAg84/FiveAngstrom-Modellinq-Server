<script setup>
const emit = defineEmits(["start-viewer"]);
const props = defineProps({
    isSidebarOpen: {
        type: Boolean,
        required: true,
    },
});

const handleStart = () => {
    emit("start-viewer");
    emit("toggle-sidebar");
};
</script>

<template>
    <div class="screen">
        <div class="container">
            <div class="home-container">
                <h1>DNCS</h1>
                <h4 class="subtitle">Digital Nets Conformational Sampling</h4>
                <button class="start-button" @click="handleStart">
                    Start <i class="fas fa-arrow-right"></i>
                </button>
            </div>
        </div>
    </div>
</template>

<style scoped>
.container {
    right: 10rem;
    align-items: center;
    display: flex;
    justify-content: center;
    padding: 10rem;
    box-sizing: border-box;
    background-color: var(--background-color);
}

@keyframes gradientShift {
    0% {
        background-position: 0% 50%;
    }
    50% {
        background-position: 100% 50%;
    }
    100% {
        background-position: 0% 50%;
    }
}

.home-container {
    width: 100%;
    height: 100%;
    text-align: center;
}

.home-container h1 {
    font-size: 3rem;
    background: linear-gradient(
        135deg,
        var(--primary-color),
        var(--secondary-color)
    );
    -webkit-background-clip: text;
    background-clip: text;
    color: transparent;
    font-weight: bold;
    text-shadow: 0px 0px 1px rgba(0, 0, 0, 0.1);
    transition: transform 0.3s ease;
}

.home-container h1:hover {
    transform: scale(1.05);
}

.subtitle {
    font-size: 1.2rem;
    color: var(--text-color);
    opacity: 0.8;
    margin-bottom: 3rem;
}

.start-button {
    display: inline-flex;
    align-items: center;
    gap: 0.5rem;
    padding: 1rem 2rem;
    font-size: 1.1rem;
    border-radius: 30px;
    background: linear-gradient(
        135deg,
        var(--primary-color),
        var(--secondary-color)
    );
    color: var(--button-text-color);
    border: none;
    cursor: pointer;
    transition:
        transform 0.3s ease,
        box-shadow 0.3s ease;
}

.start-button:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 12px var(--shadow-color);
}
</style>
