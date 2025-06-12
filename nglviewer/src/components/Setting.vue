<template>
    <div class="settings-sidebar">
        <div class="settings-header">
            <button class="close-button" @click="$emit('close-settings')">
                <i class="fas fa-times"></i>
            </button>
        </div>
        <div class="settings-content">
            <div class="representation-settings">
                <h3>Representation Settings</h3>
                <div
                    v-for="representation in representations"
                    :key="representation.name"
                    class="representation-item"
                >
                    <label class="checkbox-label">
                        <input
                            type="checkbox"
                            :value="representation.name"
                            v-model="selectedRepresentations"
                            @change="updateRepresentation(representation.name)"
                        />
                        <span class="representation-name">{{
                            representation.name
                        }}</span>
                    </label>
                </div>
            </div>
        </div>
    </div>
</template>

<script setup>
import { ref, onMounted, watch } from "vue";
const emit = defineEmits(["close-settings", "update-visibility"]);

defineProps({
    isSidebarOpen: {
        type: Boolean,
        required: true,
    },
});
const representations = ref([
    { name: "backbone", visible: false },
    { name: "cartoon", visible: false },
    { name: "line", visible: true },
    { name: "ball+stick", visible: false },
    { name: "label", visible: false },
]);

const selectedRepresentations = ref(
    representations.value
        .filter((representation) => representation.visible)
        .map((representation) => representation.name),
);

const updateRepresentation = (representationName) => {
    const isVisible =
        selectedRepresentations.value.includes(representationName);
    emit("update-visibility", representationName, isVisible);
};
</script>

<style scoped>
.settings-sidebar {
    position: fixed;
    top: 0;
    right: 0;
    height: 100vh;
    width: 20%;
    background-color: var(--sidebar-background);
    box-shadow: -2px 0 5px var(--shadow-color);
    z-index: 1000;
}

.settings-header {
    padding: 0.5rem;
    display: flex;
    justify-content: flex-end;
    border-bottom: 1px solid #ddd;
}

.settings-content {
    padding: 1rem;
}

.close-button {
    background: none;
    border: none;
    padding: 0.25rem;
    cursor: pointer;
    font-size: 1.2em;
    color: var(--text-color);
}
.representation-settings {
    padding: 1rem;
}

.representation-item {
    margin-bottom: 0.5rem;
}

.checkbox-label {
    display: flex;
    align-items: center;
    cursor: pointer;
}
.checkbox-label input[type="checkbox"] {
    margin-right: 0.5rem;
}

.representation-name {
    font-size: 1rem;
    color: var(--text-color);
}
</style>
