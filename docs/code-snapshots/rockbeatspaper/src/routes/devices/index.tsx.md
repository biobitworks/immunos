---
source: /Users/byron/projects/rockbeatspaper/src/routes/devices/index.tsx
relative: rockbeatspaper/src/routes/devices/index.tsx
generated_at: 2025-12-23 10:28
---

```tsx
import { createFileRoute } from "@tanstack/react-router";
import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import { useState } from "react";
import { useForm } from "react-hook-form";
import toast from "react-hot-toast";
import { Plus, Trash2, X } from "lucide-react";
import { ProtectedLayout } from "~/components/ProtectedLayout";
import { Card } from "~/components/Card";
import { Button } from "~/components/Button";
import { Input } from "~/components/Input";
import { useTRPC } from "~/trpc/react";
import { useAuthStore } from "~/stores/auth";

export const Route = createFileRoute("/devices/")({
  component: Devices,
});

function Devices() {
  const trpc = useTRPC();
  const queryClient = useQueryClient();
  const { authToken } = useAuthStore();
  const [isAddingDevice, setIsAddingDevice] = useState(false);

  const devicesQuery = useQuery(
    trpc.getDevices.queryOptions(
      { authToken: authToken! },
      { enabled: !!authToken }
    )
  );

  const {
    register,
    handleSubmit,
    reset,
    formState: { errors },
  } = useForm<{
    name: string;
    type: "thermostat" | "hvac" | "heater" | "cooler";
    location: string;
  }>();

  const createDeviceMutation = useMutation(
    trpc.createDevice.mutationOptions({
      onSuccess: () => {
        queryClient.invalidateQueries({ queryKey: trpc.getDevices.queryKey() });
        toast.success("Device added successfully!");
        reset();
        setIsAddingDevice(false);
      },
      onError: (error) => {
        toast.error(error.message || "Failed to add device");
      },
    })
  );

  const deleteDeviceMutation = useMutation(
    trpc.deleteDevice.mutationOptions({
      onSuccess: () => {
        queryClient.invalidateQueries({ queryKey: trpc.getDevices.queryKey() });
        toast.success("Device deleted successfully!");
      },
      onError: (error) => {
        toast.error(error.message || "Failed to delete device");
      },
    })
  );

  const onSubmit = (data: {
    name: string;
    type: "thermostat" | "hvac" | "heater" | "cooler";
    location: string;
  }) => {
    createDeviceMutation.mutate({
      authToken: authToken!,
      ...data,
    });
  };

  const handleDelete = (deviceId: number) => {
    if (confirm("Are you sure you want to delete this device?")) {
      deleteDeviceMutation.mutate({
        authToken: authToken!,
        deviceId,
      });
    }
  };

  const devices = devicesQuery.data?.devices || [];

  return (
    <ProtectedLayout>
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="flex items-center justify-between mb-8">
          <div>
            <h1 className="text-3xl font-bold text-gray-900 dark:text-white mb-2">
              Devices
            </h1>
            <p className="text-gray-600 dark:text-gray-400">
              Manage your thermal devices
            </p>
          </div>
          <Button onClick={() => setIsAddingDevice(true)}>
            <Plus className="w-4 h-4 mr-2" />
            Add Device
          </Button>
        </div>

        {/* Add Device Modal */}
        {isAddingDevice && (
          <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50 p-4">
            <Card className="w-full max-w-md">
              <div className="flex items-center justify-between mb-6">
                <h2 className="text-2xl font-bold text-gray-900 dark:text-white">
                  Add New Device
                </h2>
                <button
                  onClick={() => setIsAddingDevice(false)}
                  className="text-gray-500 hover:text-gray-700 dark:text-gray-400 dark:hover:text-gray-200"
                >
                  <X className="w-6 h-6" />
                </button>
              </div>

              <form onSubmit={handleSubmit(onSubmit)} className="space-y-4">
                <Input
                  label="Device Name"
                  {...register("name", { required: "Name is required" })}
                  error={errors.name?.message}
                  placeholder="e.g., Living Room Thermostat"
                />

                <div>
                  <label className="block text-sm font-medium text-gray-700 dark:text-gray-300 mb-1">
                    Device Type
                  </label>
                  <select
                    {...register("type", { required: "Type is required" })}
                    className="w-full px-4 py-2 border border-gray-300 dark:border-gray-600 rounded-lg focus:ring-2 focus:ring-blue-500 focus:border-transparent bg-white dark:bg-gray-700 text-gray-900 dark:text-gray-100"
                  >
                    <option value="">Select type</option>
                    <option value="thermostat">Thermostat</option>
                    <option value="hvac">HVAC System</option>
                    <option value="heater">Heater</option>
                    <option value="cooler">Cooler</option>
                  </select>
                  {errors.type && (
                    <p className="mt-1 text-sm text-red-600">
                      {errors.type.message}
                    </p>
                  )}
                </div>

                <Input
                  label="Location"
                  {...register("location", { required: "Location is required" })}
                  error={errors.location?.message}
                  placeholder="e.g., Living Room"
                />

                <div className="flex gap-3">
                  <Button
                    type="submit"
                    disabled={createDeviceMutation.isPending}
                    className="flex-1"
                  >
                    {createDeviceMutation.isPending ? "Adding..." : "Add Device"}
                  </Button>
                  <Button
                    type="button"
                    variant="outline"
                    onClick={() => setIsAddingDevice(false)}
                    className="flex-1"
                  >
                    Cancel
                  </Button>
                </div>
              </form>
            </Card>
          </div>
        )}

        {/* Devices Grid */}
        {devicesQuery.isLoading ? (
          <Card>
            <p className="text-gray-600 dark:text-gray-400">Loading devices...</p>
          </Card>
        ) : devices.length === 0 ? (
          <Card>
            <div className="text-center py-12">
              <div className="w-16 h-16 bg-gradient-to-r from-blue-600 to-purple-600 rounded-full flex items-center justify-center mx-auto mb-4">
                <Plus className="w-8 h-8 text-white" />
              </div>
              <h3 className="text-xl font-bold text-gray-900 dark:text-white mb-2">
                No devices yet
              </h3>
              <p className="text-gray-600 dark:text-gray-400 mb-6">
                Add your first device to start monitoring energy usage
              </p>
              <Button onClick={() => setIsAddingDevice(true)}>
                Add Device
              </Button>
            </div>
          </Card>
        ) : (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
            {devices.map((device) => (
              <Card key={device.id}>
                <div className="flex items-start justify-between mb-4">
                  <div>
                    <h3 className="text-lg font-bold text-gray-900 dark:text-white mb-1">
                      {device.name}
                    </h3>
                    <p className="text-sm text-gray-600 dark:text-gray-400">
                      {device.type}
                    </p>
                  </div>
                  <span
                    className={`px-3 py-1 text-xs font-semibold rounded-full ${
                      device.status === "active"
                        ? "bg-green-100 text-green-800"
                        : "bg-gray-100 text-gray-800"
                    }`}
                  >
                    {device.status}
                  </span>
                </div>

                <div className="space-y-2 mb-4">
                  <div className="flex items-center text-sm text-gray-600 dark:text-gray-400">
                    <span className="font-medium mr-2">Location:</span>
                    {device.location}
                  </div>
                  <div className="flex items-center text-sm text-gray-600 dark:text-gray-400">
                    <span className="font-medium mr-2">Added:</span>
                    {new Date(device.createdAt).toLocaleDateString()}
                  </div>
                </div>

                <Button
                  variant="danger"
                  size="sm"
                  onClick={() => handleDelete(device.id)}
                  disabled={deleteDeviceMutation.isPending}
                  className="w-full"
                >
                  <Trash2 className="w-4 h-4 mr-2" />
                  Delete
                </Button>
              </Card>
            ))}
          </div>
        )}
      </div>
    </ProtectedLayout>
  );
}

```
