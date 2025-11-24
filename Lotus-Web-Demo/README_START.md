# 如何启动服务器

## 方法 1: 使用 run.sh（推荐）

```bash
cd web_demo
./run.sh
```

## 方法 2: 手动启动

```bash
cd web_demo
source venv/bin/activate  # 如果没有 venv，先运行: python3 -m venv venv
pip install -r requirements.txt
python3 app.py
```

## 查看日志

服务器启动后，终端会显示：
- 所有 HTTP 请求
- 错误信息
- 调试信息

如果终端没有输出，可能是：
1. 服务器已经在运行（检查端口 5000）
2. 需要重启服务器以查看新日志

## 重启服务器

如果服务器已经在运行，需要先停止：

```bash
# 查找并停止进程
lsof -ti:5000 | xargs kill -9

# 然后重新启动
cd web_demo
./run.sh
```

## 访问应用

打开浏览器访问：http://localhost:5000

## 测试上传

使用示例文件：`web_demo/sample_data.h5ad`

