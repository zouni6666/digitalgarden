---
{"dg-publish":true,"permalink":"/01 研究生/01 工作/开源llm/qwen2.5/"}
---

ollama：[mp.weixin.qq.com/s/GgXy6IdaZYZVjiOb67MqcA](https://mp.weixin.qq.com/s/GgXy6IdaZYZVjiOb67MqcA)
```shell
su
sudo -i

curl -sfL https://get.gpustack.ai | sh -
pipx ensurepath
source ~/.bashrc

gpustack start --data-dir /home/data3/zhanghao/llm/qwen/ --bootstrap-password Lmz4052021!
10.133.20.72:80
cat /var/lib/gpustack/initial_admin_password	#查看密码  Lmz4052021!
```
```shell
#docker环境配置
vim /etc/systemd/system/docker.service.d/http-proxy.conf
[Service]
Environment="HTTP_PROXY=http://127.0.0.1:7890"
Environment="HTTPS_PROXY=http://127.0.0.1:7890"
Environment="NO_PROXY=localhost,127.0.0.0/8,docker-registry.somecorporation.com"

systemctl daemon-reload
systemctl restart docker
docker run -d -p 5001:5000 --name registry registry:2

sudo snap install ollama
ollama pull qwen2.5:72b
ollama cp qwen2.5:72b localhost:5001/library/qwen2.5-72b
ollama push localhost:5001/library/qwen2.5-72b --insecure

#停止自动启动
sudo systemctl list-units | grep gpustack
sudo systemctl disable gpustack.service


vim /etc/systemd/system/ollama.service
添加Environment="OLLAMA_MODELS=/path/to/ollama/models" # 记得替换
sudo systemctl daemon-reload
sudo systemctl restart ollama.service

ollama run qwen2.5:72b
退出：/bye
ps aux | grep 'ollama'
kill -9 
```



