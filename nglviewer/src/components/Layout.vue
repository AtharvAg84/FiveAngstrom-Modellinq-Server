<template>
    <div class="layout">
        <button
            v-show="!isSidebarOpen"
            class="navigation-button"
            @click="toggleSidebar"
        >
            <i class="fas fa-bars"></i>
        </button>
        <slot></slot>
    </div>
</template>

<script setup>
import { defineEmits, defineProps } from "vue";

const emit = defineEmits(["toggle-sidebar"]);

const props = defineProps({
    isSidebarOpen: {
        type: Boolean,
        required: true,
    },
});
const toggleSidebar = () => {
    emit("toggle-sidebar");
};
</script>

<style scoped>
.layout {
    display: flex;
    flex-direction: column;
    background-color: var(--background-color);
    color: var(--text-color);
}

.navigation-button {
    position: fixed;
    top: 0.5rem;
    left: 1rem;
    background: none;
    border: none;
    font-size: 1.3em;
    cursor: pointer;
    color: var(--text-color);
    padding: 0.5rem;
    z-index: 3;
}
</style>
